"""Embed downloaded OrchidRoots thumbnails with orchid-clip-v8 and DINOv2 (pre-registration Y).

What: every image in PHOTO_DIR/thumbs (grexes) and PHOTO_DIR/species_thumbs (species parents, if present) ->
L2-normalized embeddings. v8 = open_clip ViT-L-14 with musharna/orchid-clip-v8 weights (cached);
DINOv2 = facebook/dinov2-base CLS token. Images are untrusted downloads: each is opened, verified and
converted to RGB; a file that fails is listed, not silently skipped.
Output (local only): PHOTO_DIR/emb_v8.npz, PHOTO_DIR/emb_dino.npz with ids (file stem), kind (grex|species|axis|axisgrex; axis = held-out photos of the 9 mini-catt axis species in PHOTO_DIR/axis_heldout, axisgrex = second photos of the trajectory-axis grexes in PHOTO_DIR/axis_grex), emb.
Usage: python -I analysis/embed_thumbs.py PHOTO_DIR [--limit N] [--device cuda|cpu]
"""

import glob
import json
import os
import sys

import numpy as np
import torch
from PIL import Image


def images(photo_dir, limit):
    out = []
    for kind, sub in (
        ("grex", "thumbs"),
        ("species", "species_thumbs"),
        ("axis", "axis_heldout"),
        ("axisgrex", "axis_grex"),
    ):
        for f in sorted(glob.glob(os.path.join(photo_dir, sub, "*"))):
            if f.endswith(
                ".jsonl"
            ):  # axis_heldout keeps its fetch manifest beside the photos
                continue
            out.append((kind, os.path.splitext(os.path.basename(f))[0], f))
    return out[:limit] if limit else out


def load(f):
    with Image.open(f) as im:
        im.verify()
    return Image.open(f).convert("RGB")


def main(photo_dir, limit=0, device="cuda"):
    items = images(photo_dir, limit)
    bad, ok = [], []
    for kind, i, f in items:
        try:
            ok.append((kind, i, load(f)))
        except Exception as e:  # noqa: BLE001 - reported below, run fails if any
            bad.append((f, repr(e)))
    if bad:
        raise SystemExit(f"{len(bad)} unreadable images: {bad[:5]}")
    import open_clip
    from huggingface_hub import hf_hub_download
    from transformers import AutoImageProcessor, AutoModel

    cfg = json.load(
        open(hf_hub_download("musharna/orchid-clip-v8", "model_config.json"))
    )
    model, _, prep = open_clip.create_model_and_transforms(
        cfg["model_name"], pretrained=None
    )
    state = torch.load(
        hf_hub_download("musharna/orchid-clip-v8", "open_clip_pytorch_model.bin"),
        map_location="cpu",
        weights_only=False,
    )
    model.load_state_dict(state["state_dict"])
    model.to(device).eval()
    dp = AutoImageProcessor.from_pretrained("facebook/dinov2-base")
    dm = AutoModel.from_pretrained("facebook/dinov2-base").to(device).eval()
    ev, ed = [], []
    with torch.no_grad():
        for s in range(0, len(ok), 64):
            ims = [im for _, _, im in ok[s : s + 64]]
            x = torch.stack([prep(im) for im in ims]).to(device)
            ev.append(
                torch.nn.functional.normalize(model.encode_image(x).float(), dim=-1)
                .cpu()
                .numpy()
            )
            d = (
                dm(**dp(images=ims, return_tensors="pt").to(device))
                .last_hidden_state[:, 0]
                .float()
            )
            ed.append(torch.nn.functional.normalize(d, dim=-1).cpu().numpy())
    ids = np.array([i for _, i, _ in ok])
    kind = np.array([k for k, _, _ in ok])
    for name, e in (("v8", ev), ("dino", ed)):
        e = np.concatenate(e)
        assert e.shape[0] == len(ok) and np.allclose(
            np.linalg.norm(e, axis=1), 1, atol=1e-4
        )
        np.savez(
            os.path.join(photo_dir, f"emb_{name}{'_TEST' if limit else ''}.npz"),
            ids=ids,
            kind=kind,
            emb=e,
        )
        print(
            f"{name}: {e.shape} ({(kind == 'grex').sum()} grex, {(kind == 'species').sum()} species)"
        )


if __name__ == "__main__":
    a = sys.argv[1:]
    lim = int(a[a.index("--limit") + 1]) if "--limit" in a else 0
    dev = a[a.index("--device") + 1] if "--device" in a else "cuda"
    main(a[0], lim, dev)
