"""Extract the figure bitmaps of arXiv:2609.12184 that contain the measured transfer curves.

The paper PDF is not redistributed in this repository. Download it from
https://arxiv.org/abs/2609.12184 (v2), save it as data/refs/jeong2026_agentic_tcad.pdf and run

    python scripts/extract_figs.py
    python scripts/digitize_paper.py

This writes data/figs/p5_img0.png (Fig. 5), p5_img1.png (Fig. 6) and p6_img0.png (Fig. 7), the
inputs of digitize_paper.py. The digitised result (data/meas.csv) is already included, so this
step is only needed to re-run or audit the digitisation.
"""
from __future__ import annotations

import argparse
import hashlib
import os

# MD5 of the bitmaps extracted from the PDF version used for data/meas.csv (arXiv v2)
EXPECTED_MD5 = {
    "p5_img0.png": "a59e3b354b8b4f55b2899fec152924d3",
    "p5_img1.png": "ec64f99da2a4a4f9ae9143cd172dce94",
    "p6_img0.png": "76390d5bc70f662271671aab1afde7f9",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default="data/refs/jeong2026_agentic_tcad.pdf")
    ap.add_argument("--out", default="data/figs")
    args = ap.parse_args()
    import pymupdf  # pip install pymupdf
    os.makedirs(args.out, exist_ok=True)
    doc = pymupdf.open(args.pdf)
    for pno in (4, 5):                              # 0-based: pages 5 and 6 hold Figs. 5-7
        page = doc[pno]
        for k, img in enumerate(page.get_images(full=True)):
            info = doc.extract_image(img[0])
            fn = os.path.join(args.out, f"p{pno + 1}_img{k}.{info['ext']}")
            with open(fn, "wb") as f:
                f.write(info["image"])
            md5 = hashlib.md5(info["image"]).hexdigest()
            want = EXPECTED_MD5.get(os.path.basename(fn))
            check = "" if want is None else ("  [MD5 OK]" if md5 == want else "  [MD5 differs: not the PDF version used here]")
            print(f"{fn}: {info['width']} x {info['height']}{check}")


if __name__ == "__main__":
    main()
