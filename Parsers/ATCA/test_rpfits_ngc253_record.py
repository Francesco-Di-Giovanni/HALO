import ctypes as C

LIB_PATH = "/home/franz/miniconda3/envs/rpfits_env/lib/librpfits.so"
FITS_PATH = "/home/franz/FRO/Parsers/ATCA/Surveys/2002-02-09_2302.C1025"

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


class IfT(C.Structure):
    _fields_ = [
        ("n_if", C.c_int),
        ("if_invert", C.c_int * MAX_IF),
        ("if_nfreq", C.c_int * MAX_IF),
        ("if_nstok", C.c_int * MAX_IF),
        ("if_sampl", C.c_int * MAX_IF),
        ("if_found", C.c_int),
        ("if_num", C.c_int * MAX_IF),
        ("if_simul", C.c_int * MAX_IF),
        ("if_chain", C.c_int * MAX_IF),
    ]


class DoublesT(C.Structure):
    _fields_ = [
        ("axis_offset", C.c_double * ANT_MAX),
        ("dec", C.c_double),
        ("dfreq", C.c_double),
        ("cu_cal1", C.c_double * 32),
        ("cu_cal2", C.c_double * 32),
        ("cu_ut", C.c_double * 32),
        ("feed_cal", C.c_double * (ANT_MAX * MAX_IF * 8)),
        ("feed_pa", C.c_double * (2 * ANT_MAX)),
        ("fg_ut", C.c_double * (2 * MAX_FG)),
        ("freq", C.c_double),
        ("if_bw", C.c_double * MAX_IF),
        ("if_ref", C.c_double * MAX_IF),
        ("if_freq", C.c_double * MAX_IF),
        ("mt_humid", C.c_double * 256),
        ("mt_press", C.c_double * 256),
        ("mt_temp", C.c_double * 256),
        ("mt_ut", C.c_double * 256),
        ("nx_ut", C.c_double * MAX_NX),
        ("ra", C.c_double),
        ("rfreq", C.c_double),
        ("rp_c", C.c_double * 12),
        ("rp_djmrefp", C.c_double),
        ("rp_djmreft", C.c_double),
        ("rp_utcmtai", C.c_double),
        ("su_dec", C.c_double * MAX_SU),
        ("su_ra", C.c_double * MAX_SU),
        ("su_rad", C.c_double * MAX_SU),
        ("su_decd", C.c_double * MAX_SU),
        ("su_pra", C.c_double * MAX_SU),
        ("su_pdec", C.c_double * MAX_SU),
        ("su_prad", C.c_double * MAX_SU),
        ("su_pdecd", C.c_double * MAX_SU),
        ("vel1", C.c_double),
        ("x", C.c_double * ANT_MAX),
        ("x_array", C.c_double),
        ("y", C.c_double * ANT_MAX),
        ("y_array", C.c_double),
        ("z", C.c_double * ANT_MAX),
        ("z_array", C.c_double),
    ]


param_ = ParamT.in_dll(lib, "param_")
names_ = NamesT.in_dll(lib, "names_")
if_ = IfT.in_dll(lib, "if_")
doubles_ = DoublesT.in_dll(lib, "doubles_")

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

# Skip through headers until we reach NGC253 (scan 4)
target_scan = 4
scan_num = 0
while scan_num < target_scan:
    r = call_rpfitsin(-1)
    scan_num += 1
    obj = names_.object.decode(errors="replace").strip()
    print(f"Header {scan_num}: object={obj!r} jstat={r}")
    if obj == "NGC253":
        break
    # Skip this scan's data
    while True:
        r2 = call_rpfitsin(0)
        if r2 in (1, 3):
            break

print("=== NGC253 header reached ===")
print("n_if:", if_.n_if)
for i in range(if_.n_if):
    print(f"  IF {i}: if_num={if_.if_num[i]} if_nfreq={if_.if_nfreq[i]} "
          f"if_nstok={if_.if_nstok[i]} if_freq={doubles_.if_freq[i]/1e6:.4f} MHz "
          f"if_bw={doubles_.if_bw[i]/1e6:.4f} MHz if_invert={if_.if_invert[i]} "
          f"if_simul={if_.if_simul[i]}")
print("RA (rad):", doubles_.ra, " Dec (rad):", doubles_.dec)
print("nstok:", param_.nstok, "nfreq:", param_.nfreq, "intime:", param_.intime)
print("bunit:", names_.bunit)
print("datobs:", names_.datobs)

print()
print("=== First data record of NGC253 ===")
r2 = call_rpfitsin(0)
print("jstat:", r2)
print("baseline:", baseline.value)
print("if_no:", if_no.value)
print("sourceno:", sourceno.value)
print("flag:", flag.value)
print("bin:", bin_.value)
print("ut (s):", ut.value)
print("u,v,w:", u.value, v.value, w.value)
print("First 10 vis values:", list(vis[:10]))
print("First 10 wgt values:", list(wgt[:10]))

print()
print("=== Next 15 records (baseline/if_no/flag/sourceno/ut) ===")
for i in range(15):
    r3 = call_rpfitsin(0)
    if r3 != 0:
        print(f"  record {i}: jstat={r3} (stop)")
        break
    print(f"  record {i}: baseline={baseline.value} if_no={if_no.value} flag={flag.value} sourceno={sourceno.value} ut={ut.value:.3f}")

call_rpfitsin(1)
print("Closed file.")
