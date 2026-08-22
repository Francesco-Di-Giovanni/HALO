import ctypes as C
import csv

LIB_PATH = "/home/franz/miniconda3/envs/rpfits_env/lib/librpfits.so"
FITS_PATH = "/home/franz/FRO/Parsers/ATCA/Surveys/2002-02-09_2302.C1025"
OUT_CSV = "/home/franz/FRO/Parsers/ATCA/Surveys/ngc253_baseline258_spectrum.csv"
OUT_PNG = "/home/franz/FRO/Parsers/ATCA/Surveys/ngc253_baseline258_spectrum.png"

ANT_MAX = 16
MAX_IF = 48
MAX_SU = 2048
MAX_CARD = 2240
MAX_FG = 32
MAX_NX = 256

lib = C.CDLL(LIB_PATH)


class ParamT(C.Structure):
    _fields_ = [
        ("nstok", C.c_int), ("nfreq", C.c_int), ("ncount", C.c_int),
        ("intime", C.c_int), ("nscan", C.c_int), ("write_wt", C.c_int),
        ("ncard", C.c_int), ("intbase", C.c_float), ("data_format", C.c_int),
    ]


class NamesT(C.Structure):
    _fields_ = [
        ("object", C.c_char * 16), ("instrument", C.c_char * 16),
        ("cal", C.c_char * 16), ("rp_observer", C.c_char * 16),
        ("datobs", C.c_char * 12), ("datwrit", C.c_char * 12),
        ("file", C.c_char * 256), ("datsys", C.c_char * 8),
        ("version", C.c_char * 20), ("coord", C.c_char * 8),
        ("sta", C.c_char * (ANT_MAX * 8)),
        ("feed_type", C.c_char * (2 * ANT_MAX * 2)),
        ("card", C.c_char * (MAX_CARD * 80)),
        ("if_cstok", C.c_char * (4 * MAX_IF * 2)),
        ("su_name", C.c_char * (MAX_SU * 16)),
        ("su_cal", C.c_char * (MAX_SU * 4)),
        ("fg_reason", C.c_char * (MAX_FG * 24)),
        ("nx_source", C.c_char * (MAX_NX * 16)),
        ("nx_date", C.c_char * (MAX_NX * 12)),
        ("rpfitsversion", C.c_char * 20), ("bunit", C.c_char * 16),
        ("obstype", C.c_char * 16), ("errmsg", C.c_char * 80),
    ]


param_ = ParamT.in_dll(lib, "param_")
names_ = NamesT.in_dll(lib, "names_")

names_.file = FITS_PATH.encode("ascii")

MAXBUF = 500000
vis = (C.c_float * MAXBUF)()
wgt = (C.c_float * MAXBUF)()
baseline = C.c_int(0)
ut = C.c_float(0)
u = C.c_float(0)
v = C.c_float(0)
w = C.c_float(0)
flag = C.c_int(0)
bin_ = C.c_int(0)
if_no = C.c_int(0)
sourceno = C.c_int(0)


def call_rpfitsin(jstat_val):
    jstat = C.c_int(jstat_val)
    lib.rpfitsin_(
        C.byref(jstat), vis, wgt, C.byref(baseline), C.byref(ut),
        C.byref(u), C.byref(v), C.byref(w), C.byref(flag), C.byref(bin_),
        C.byref(if_no), C.byref(sourceno)
    )
    return jstat.value


r = call_rpfitsin(-3)
print(f"open: jstat={r}")

# Skip headers until NGC253
scan_num = 0
while True:
    r = call_rpfitsin(-1)
    scan_num += 1
    obj = names_.object.decode(errors="replace").strip()
    print(f"Header {scan_num}: object={obj!r}")
    if obj == "NGC253":
        break

print(f"Reached NGC253. nstok={param_.nstok} nfreq={param_.nfreq}")

TARGET_BASELINE = 258
found = False
while True:
    r2 = call_rpfitsin(0)
    if r2 == 1:
        print("End of scan reached without finding target baseline.")
        break
    if r2 == 3:
        print("EOF reached.")
        break
    if r2 != 0:
        continue  # intermediate table, keep going
    if baseline.value == TARGET_BASELINE and if_no.value == 1:
        found = True
        break

if not found:
    print("Target baseline not found in this scan.")
    call_rpfitsin(1)
    raise SystemExit

nstok = param_.nstok
nfreq = param_.nfreq
print(f"Found baseline={baseline.value} if_no={if_no.value} ut={ut.value:.3f} bin={bin_.value} flag={flag.value}")

# Extract complex vis: assume ordering index = 2*(ch*nstok + pol) + {0,1}
amp = [[0.0] * nfreq for _ in range(nstok)]
phase = [[0.0] * nfreq for _ in range(nstok)]
weight_out = [[0.0] * nfreq for _ in range(nstok)]

import math
for ch in range(nfreq):
    for pol in range(nstok):
        idx = ch * nstok + pol
        re = vis[2 * idx]
        im = vis[2 * idx + 1]
        amp[pol][ch] = math.sqrt(re * re + im * im)
        phase[pol][ch] = math.degrees(math.atan2(im, re))
        weight_out[pol][ch] = wgt[idx]

# Save CSV
with open(OUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    header = ["channel"]
    for pol in range(nstok):
        header += [f"amp_pol{pol}", f"phase_pol{pol}_deg", f"weight_pol{pol}"]
    writer.writerow(header)
    for ch in range(nfreq):
        row = [ch]
        for pol in range(nstok):
            row += [amp[pol][ch], phase[pol][ch], weight_out[pol][ch]]
        writer.writerow(row)

print(f"CSV saved: {OUT_CSV}")

# Plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
for pol in range(nstok):
    axes[0].plot(range(nfreq), amp[pol], label=f"Pol {pol}")
    axes[1].plot(range(nfreq), phase[pol], label=f"Pol {pol}")
axes[0].set_ylabel("Amplitude (Jy, uncalibrated)")
axes[0].set_title(f"NGC253 - baseline {TARGET_BASELINE} (ant 1-2) - IF0 1419 MHz - raw visibility")
axes[0].legend()
axes[1].set_ylabel("Phase (deg)")
axes[1].set_xlabel("Channel")
axes[1].legend()
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=120)
print(f"PNG saved: {OUT_PNG}")

call_rpfitsin(1)
print("Closed file.")
