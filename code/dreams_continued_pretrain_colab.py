"""dreams_continued_pretrain_colab.py — GPU step (run on Colab, NOT this Mac).

Continued masked-peak self-supervision of DreaMS on our domain corpus, to make the model
speak Orbitrap-HILIC and (hopefully) pull the rich-but-isolated uncurated bins into structure.

HOW TO RUN (Colab, GPU runtime = T4 is enough):
  1. New Colab notebook → Runtime → Change runtime type → GPU.
  2. Upload two files (Files pane) to /content/:
        data/dreams_ssl/corpus.hdf5         (15,051 spectra, the SSL corpus)
        data/dreams_uncurated/probe.mgf     (the 4,828-spectrum probe, for the re-test)
  3. Upload this script and run:  !python dreams_continued_pretrain_colab.py
     (or paste the STEP blocks into cells one at a time — easier to debug.)
  4. Download the two outputs it writes:
        emb_base_backbone.npy      (baseline = pretrained SSL backbone)
        emb_adapted_backbone.npy   (after domain adaptation)
     Bring them back to data/dreams_uncurated/ and run the local Phase-3 compare.

NOTES / KNOWN RISKS (first run is experimental — paste any error back and we iterate):
  • Resume patch: train.py only wired pre_trained_pth for fine-tuning; we patch the
    pre-training branch to DreaMS.load_from_checkpoint(...). Architecture args MUST match the
    released ssl_model.ckpt → we copy pre_train.sh's arch flags verbatim.
  • Small corpus (15k vs 700M) → LOW lr (1e-5) + FEW epochs (15) to adapt, not overwrite.
  • ret_order_loss_w=0 (our hdf5 has no RT-pairing groups) → pure masked-peak objective.
  • Embeddings for the compare come from the BACKBONE via dreams_intermediates(precursor_only),
    NOT dreams_embeddings (that's the ContrastiveHead = a different, fine-tuned model).
    So the baseline is the ORIGINAL backbone, isolating the adaptation effect.
"""
import os, re, subprocess, sys, glob
from pathlib import Path

CORPUS = '/content/corpus.hdf5'
PROBE  = '/content/probe.mgf'
RUN    = 'domain_adapt'                       # project/job key
ARCH = [  # verbatim from DreaMS/dreams/training/pre_train.sh — must match ssl_model.ckpt
    '--dformat', 'A', '--model', 'DreaMS',
    '--n_layers', '7', '--n_heads', '8', '--d_peak', '44', '--d_fourier', '980',
    '--ff_peak_depth', '1', '--ff_fourier_depth', '5', '--ff_fourier_d', '512', '--ff_out_depth', '1',
    '--prec_intens', '1.1', '--max_peaks_n', '60', '--attn_mech', 'dot-product',
    '--fourier_strategy', 'lin_float_int', '--focal_loss_gamma', '5',
    '--no_transformer_bias', '--pre_norm', '--graphormer_mz_diffs',
    '--att_dropout', '0.1', '--residual_dropout', '0.1', '--ff_dropout', '0.1', '--dropout', '0.1',
    '--train_objective', 'mask_mz_hot', '--hot_mz_bin_size', '0.05',
    '--frac_masks', '0.3', '--mask_peaks', '--mask_intens_strategy', 'intens_p',
]


def sh(cmd):
    print('+', cmd); subprocess.run(cmd, shell=True, check=True)


# ── STEP 1 — install + clone ───────────────────────────────────────────────
def setup():
    sh('pip -q install uv || true')
    if not Path('/content/DreaMS').exists():
        sh('git clone --depth 1 https://github.com/pluskal-lab/DreaMS.git /content/DreaMS')
        sh('pip -q install -e /content/DreaMS')
        sh('pip -q install ms_entropy h5py')
    # PATCH 1: enable continued pretraining + strict=False (ssl_model.ckpt carries a
    #          retention-order head ro_out we don't build when ret_order_loss_w=0).
    tp = Path('/content/DreaMS/dreams/training/train.py')
    s = tp.read_text()
    if 'PATCH: continued pretraining' not in s:
        s = s.replace(
            "        if args.model == 'DreaMS':\n"
            "            if not args.pre_trained_pth:\n"
            "            #     model = DreaMS.load_from_checkpoint(args.pre_trained_pth, map_location=torch.device(device))\n"
            "            # else:\n"
            "                model = DreaMS(args, spec_preproc)",
            "        if args.model == 'DreaMS':\n"
            "            if args.pre_trained_pth:   # PATCH: continued pretraining\n"
            "                model = DreaMS.load_from_checkpoint(\n"
            "                    args.pre_trained_pth, map_location=torch.device(device),\n"
            "                    args=args, spec_preproc=spec_preproc, strict=False)\n"
            "            else:\n"
            "                model = DreaMS(args, spec_preproc)")
        tp.write_text(s)
        assert 'PATCH: continued pretraining' in s, 'PATCH 1 FAILED — patch train.py block by hand'
        print('patched train.py (resume + strict=False)')

    # PATCH 2: upstream bug — to_classes() has its return_num_classes 2-tuple return commented
    #          out, but to_hot() unpacks it (breaks the mask_mz_hot objective).
    sp = Path('/content/DreaMS/dreams/utils/spectra.py')
    t = sp.read_text()
    if '#     return classes, num_classes + len(special_vals)' in t:
        t = t.replace(
            "    # if return_num_classes:\n"
            "    #     return classes, num_classes + len(special_vals)\n"
            "    return classes",
            "    if return_num_classes:   # PATCH: to_hot relies on this 2-tuple return\n"
            "        return classes, num_classes + len(special_vals)\n"
            "    return classes")
        sp.write_text(t)
        print('patched spectra.py (to_classes return_num_classes)')


# ── STEP 2 — fetch the pretrained SSL backbone ─────────────────────────────
def get_ckpt():
    sys.path.insert(0, '/content/DreaMS')
    from dreams.definitions import PRETRAINED
    from dreams import utils
    ck = Path(PRETRAINED) / 'ssl_model.ckpt'
    if not ck.exists():
        utils.io.download_pretrained_model('ssl_model.ckpt')  # may be utils.download_pretrained_model
    print('ssl backbone:', ck, ck.exists())
    return str(ck)


# ── STEP 3 — continued pretraining ─────────────────────────────────────────
def train(ssl_ckpt):
    cmd = [sys.executable, '/content/DreaMS/dreams/training/train.py',
           '--project_name', RUN, '--job_key', RUN, '--run_name', RUN,
           '--train_regime', 'pre-training', '--pre_trained_pth', ssl_ckpt,
           '--dataset_pth', CORPUS,
           '--lr', '1e-5', '--max_epochs', '15', '--batch_size', '32',
           '--num_devices', '1', '--ret_order_loss_w', '0',
           '--val_check_interval', '1.0', '--n_warmup_steps', '200',
           '--train_precision', '32', '--no_wandb', '--save_top_k', '-1', '--seed', '42'] + ARCH
    print('+', ' '.join(cmd)); subprocess.run(cmd, check=True)
    cks = sorted(glob.glob(f'/content/{RUN}/**/*.ckpt', recursive=True), key=os.path.getmtime)
    assert cks, 'no checkpoint produced — check training logs'
    print('adapted ckpt:', cks[-1]); return cks[-1]


# ── STEP 4 — backbone embeddings: baseline vs adapted (apples-to-apples) ────
def embed(ssl_ckpt, adapted_ckpt):
    import numpy as np
    from dreams.api import dreams_intermediates
    for tag, ck in [('base', ssl_ckpt), ('adapted', adapted_ckpt)]:
        emb = np.asarray(dreams_intermediates(ck, PROBE, precursor_only=True, batch_size=32))
        out = f'/content/emb_{tag}_backbone.npy'; np.save(out, emb)
        print(f'wrote {out}  shape={emb.shape}')


if __name__ == '__main__':
    setup()
    ssl_ckpt = get_ckpt()
    adapted = train(ssl_ckpt)
    embed(ssl_ckpt, adapted)
    print('\nDONE — download emb_base_backbone.npy + emb_adapted_backbone.npy to data/dreams_uncurated/')
