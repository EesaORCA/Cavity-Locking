#%%
"""
Red Pitaya cavity lock  (single PID, lock-in / PDH-style error signal)
Out1 -> Laser AC modulation input (could be swapped to an EOM)
Out2 -> PDm200 amplifier -> piezo
PDA10CS2 Photodiode 0dB gain -> In2 

Iq module used for fast 10MHz modulation and demod of the laser signal. 
PID module used both for locking and for scanning in mode where p=i=0
Scope module used to visualise traces and monitor transmission during the lock.




"""

import subprocess
import time
import logging
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy.signal import savgol_filter, find_peaks

from pyrpl import Pyrpl
# ===========================  PARAMETERS  ============================
# ---- connection -----------------------------------------------------
HOSTNAME        = "rp-f0efad.local"
CONFIG          = "PIDlock"
RELOAD_FPGA     = True        #SSH Reload
FPGA_ATTEMPTS   = 3

# ---- modulation / demodulation (iq0) --------------------------------
F_MOD           = 10e6        # Hz, modulation frequency
MOD_DEPTH       = 0.0032       # V, modulation amplitude on OUT1

PHASE           = 120         # deg, demodulation phase
IQ_BANDWIDTH    = 2e3         # Hz, low-pass on demodulated output 

AC_BANDWIDTH    = 2e4         # Hz, ac-coupling high-pass on the iq input
QUAD_FACTOR     = 1           # Digital amplification
IQ_INPUT        = 'in2'       # cavity transmission photodiode
IQ_OUTPUT       = 'out1'      # where the modulation tone goes

# ---- actuator output  -----------

V_MIN           = 0.0         # V, hard FPGA clamp, low, ensures no negative voltage is sent to piezo.
V_MAX           = 1.0         # V, hard FPGA clamp, high
V_GUESS         = 0.7       # V, expected resonance position

# ---- acquisition sweep ----------------------------------------------
ACQ_SPAN        = 0.1       # V, +/- span for the FIRST search
RELOCK_SPAN     = 0.02        # V, +/- span when re-acquiring
SWEEP_STEP      = 2.5e-4        # V, step size
SWEEP_SETTLE    = 3e-3        # s, dwell after each ival write
SMOOTH_POINTS   = 5           # savgol window for peak finding (odd >=5, or 0)
FIT_POINTS      = 2           # points either side of peak for the slope fit

# ---- resonance acceptance criteria ----------------------------------
MIN_PEAK        = 0.2         # V, minimum peak height for locking peak
MAX_CANDIDATES  = 2           # how many peaks to try, best transmission first
PEAK_SEPARATION = 0.002       # V, minimum spacing between distinct candidates

MIN_ERR_AMP     = 0.03        # V, minimum differential size

ERR_WINDOW      = 0.005       # V, range searched around a resonance
REQUIRE_BRACKET = True        # extrema must straddle the transmission peak
MIN_SLOPE       = 20.0        # V/V, minimum |discriminant slope| at the peak

# ---- servo ----------------------------------------------------------
I_GAIN          = 20        # Hz, integrator unity-gain frequency (magnitude)
                              #   f_unity ~= I_GAIN * |slope|
P_GAIN          = 0.0         # proportional gain (magnitude)
SIGN_OVERRIDE   = -1          # None = from measured slope, or force +1 / -1
SETPOINT_FROM_SCAN = True     # lock where the error sat at max transmission

# ---- network analyser (plant measurement while locked) --------------
RUN_NA_AT_START = False      # measure the plant right after the lock settles
NA_START        = 20.0        # Hz
NA_STOP         = 2e4       # Hz
NA_POINTS       = 201
NA_RBW          = 100       # Hz, resolution bandwidth (lower = slower, cleaner)
NA_AVG          = 1           # averages per point
NA_AMPLITUDE    = 3e-4        # V at OUT2. MUST stay well inside the linear
                              #   region: half-linewidth is ~1.2 mV. Halve it
                              #   and check the transfer function is unchanged.
NA_I_GAIN       = 0.1         # Hz, reduced integrator gain during the sweep so
                              #   the loop is transparent above ~20 Hz but
                              #   still holds against drift
NA_IQ_BANDWIDTH = 20e3        # Hz, widen the demodulator so its own poles are
                              #   not what you measure. Restored afterwards.
NA_SETTLE       = 3.0         # s, let the loop re-settle after changing gain
NA_LOGSCALE     = True
NA_SEARCH_MIN   = 100.0       # Hz, ignore below this when hunting a resonance
NA_TARGET_MARGIN = 4.0        # put unity gain at f_resonance / this
NA_DEEMBED      = True        # correct M -> P using the known PID transfer
NA_SAVE         = True

# ---- monitoring / watchdog ------------------------------------------
MEAS_DURATION   = 1e-3        # s, scope trace length per readback
MONITOR_DT      = 0.2         # How often scope values are monitored 
LOCK_THRESHOLD  = 0.5         # relock if transmission < this fraction of peak
RUN_TIME        = 30000        # s, total time to hold the lock
MAX_RELOCKS     = 20
PRINT_EVERY     = 1

# ---- plotting / analysis window -------------------------------------
PLOT_SETTLE     = 10.0        # s, discard this much after every engage before
                              #   plotting or computing statistics (piezo creep
                              #   dominates the first ~30 s after a sweep)

# ---- shutdown -------------------------------------------------------
RAMP_DOWN       = True
RAMP_DOWN_TIME  = 1.0
RAMP_DOWN_STEPS = 200

# ---- logging / plots ------------------------------------------------
PLOT_ACQUISITION = True
PLOT_LOG         = True
SAVE_LOG         = True
SAVE_SWEEP       = True       # save the accepted acquisition sweep too
SAVE_DIR         = r'C:\Users\Lab User 10\OneDrive - ORCA Computing\Documents\redpitayacavitydatascans'

# ---- Cleanup and initialisation-----------
try:
    _stale = p
except NameError:
    _stale = None
if _stale is not None:
    _log = logging.getLogger('pyrpl.redpitaya')
    _lvl = _log.level
    _log.setLevel(logging.CRITICAL)
    for _fn in (lambda: _stale.rp.end_all(),   
                lambda: _stale.close()):       
        try:                                   
            _fn()
        except Exception:
            pass
    _log.setLevel(_lvl)
    del _stale, p

if RELOAD_FPGA:
    for attempt in range(FPGA_ATTEMPTS):
        r = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
             f"root@{HOSTNAME}",
             "killall -q overlay.sh monitor-server; sleep 1; "
             "/opt/redpitaya/sbin/overlay.sh pyrpl >/dev/null 2>&1; "
             "cat /sys/class/fpga_manager/fpga0/state"],
            capture_output=True, text=True, timeout=45)
        if "operating" in r.stdout:
            break
        print(f"attempt {attempt+1} failed: {r.stdout.strip()} {r.stderr.strip()}")
    else:
        raise RuntimeError("FPGA reload failed after 3 attempts")


# ---- connect --------------------------------------------------------
p = Pyrpl(config=CONFIG, hostname=HOSTNAME, gui=False, reloadfpga=False)
s = p.rp.scope
iq = p.rp.iq0
asg = p.rp.asg1
pid = p.rp.pid0
sampl=p.rp.sampler


# ---- modulation -----------------------------------------------------
iq.setup(frequency=F_MOD,
         amplitude=MOD_DEPTH,
         phase=PHASE,
         bandwidth=[IQ_BANDWIDTH, IQ_BANDWIDTH],
         acbandwidth=AC_BANDWIDTH,
         gain=0,
         input=IQ_INPUT,
         output_direct=IQ_OUTPUT,
         output_signal='quadrature',
         quadrature_factor=QUAD_FACTOR)
print("iq0 bandwidth readback:", iq.bandwidth)   # FPGA quantises the corner


# ---- make sure nothing but pid0 can reach OUT2 -----------------------
asg.amplitude = 0
asg.output_direct = 'off'
for name in ('asg0', 'asg1', 'pid1', 'pid2', 'iq1', 'iq2'):
    try:
        mod = getattr(p.rp, name)
        if mod.output_direct == 'out2':
            mod.output_direct = 'off'
            print(f"  disconnected {name} from out2")
    except Exception as e:
        print(f"  {name} check warning:", e)


# ---- PID ------------------------------------------------------------
pid.input = 'iq0'
pid.output_direct = 'out2'
pid.setpoint =0
pid.p =0
pid.i =0
pid.d=0
pid.inputfilter = [0, 0, 0, 0]
pid.min_voltage = V_MIN
pid.max_voltage = V_MAX
pid.ival = float(np.clip(V_GUESS, V_MIN, V_MAX))


# ---- scope  -----------------------------------
s.setup(input1=IQ_INPUT, input2='iq0',
        duration=MEAS_DURATION,
        trigger_source='immediately',
        trace_average=1,
        rolling_mode=False,
        average=True)


# Functions

def measure():
    """Return (mean transmission, mean raw iq0 signal) in volts."""
    d = np.array(s.single())
    return float(np.mean(d[0])), float(np.mean(d[1]))


def set_output(v):
    """Hold OUT2 at v volts (only valid while p = i = 0)."""
    pid.ival = float(np.clip(v, V_MIN, V_MAX))


def sweep(center, span, step=SWEEP_STEP):
    """Open-loop sweep of OUT2. Returns (v, transmission, error)."""
    pid.p = 0.0
    pid.i = 0.0
    v = np.arange(center - span, center + span + 0.5 * step, step)
    v = v[(v >= V_MIN) & (v <= V_MAX)]
    if len(v) < 5:
        raise RuntimeError("sweep range is empty - check V_MIN/V_MAX/V_GUESS")

    trans = np.zeros(len(v))
    err = np.zeros(len(v))
    t0 = time.perf_counter()
    for n, vv in enumerate(v):
        set_output(vv)
        time.sleep(SWEEP_SETTLE)
        trans[n], err[n] = measure()
    print(f"  swept {len(v)} points over "
          f"{v[0]:.4f} - {v[-1]:.4f} V in {time.perf_counter()-t0:.2f} s")
    return v, trans, err


def smooth(trans):
    if SMOOTH_POINTS >= 5 and len(trans) > SMOOTH_POINTS:
        return savgol_filter(trans, SMOOTH_POINTS | 1, 2)
    return trans


def candidates(v, trans_s):
    """Indices of candidate peaks, best transmission first."""
    dist = max(int(round(PEAK_SEPARATION / SWEEP_STEP)), 1)
    idx, _ = find_peaks(trans_s, height=MIN_PEAK, distance=dist)
    if len(idx) == 0:                        # fall back to the global maximum
        idx = np.array([int(np.argmax(trans_s))])
    order = np.argsort(trans_s[idx])[::-1]
    return idx[order][:MAX_CANDIDATES]


def characterise(v, trans_s, err, k):
    """Measure the discriminant of the resonance at index k."""
    lo = max(k - FIT_POINTS, 0)
    hi = min(k + FIT_POINTS + 1, len(v))
    slope = float(np.polyfit(v[lo:hi], err[lo:hi], 1)[0])   # V_err per V_out

    win = np.abs(v - v[k]) <= ERR_WINDOW
    e_win, v_win = err[win], v[win]
    i_hi, i_lo = int(np.argmax(e_win)), int(np.argmin(e_win))
    e_c = float(err[k])

    d = dict(k=int(k), v_peak=float(v[k]), t_peak=float(trans_s[k]),
             e_peak=e_c, slope=slope, n_win=int(win.sum()),
             e_hi=float(e_win[i_hi]), v_hi=float(v_win[i_hi]),
             e_lo=float(e_win[i_lo]), v_lo=float(v_win[i_lo]))
    d['amp_up'] = d['e_hi'] - e_c            # excursion above the lock point
    d['amp_dn'] = e_c - d['e_lo']            # excursion below the lock point
    d['bracketed'] = (min(d['v_lo'], d['v_hi']) < d['v_peak']
                      < max(d['v_lo'], d['v_hi']))
    return d


def validate(d):
    """Return a list of reasons this resonance is not lockable ([] = good)."""
    bad = []
    if d['t_peak'] < MIN_PEAK:
        bad.append(f"transmission {d['t_peak']:.3f} V < MIN_PEAK {MIN_PEAK} V")
    if d['n_win'] < 5:
        bad.append(f"only {d['n_win']} sweep points within "
                   f"+/-{ERR_WINDOW*1e3:.1f} mV - reduce SWEEP_STEP "
                   f"or widen ERR_WINDOW")
    if MIN_ERR_AMP > 0:
        if d['amp_up'] < MIN_ERR_AMP:
            bad.append(f"error rises only {d['amp_up']*1e3:+.1f} mV above the "
                       f"lock point (need {MIN_ERR_AMP*1e3:.0f} mV)")
        if d['amp_dn'] < MIN_ERR_AMP:
            bad.append(f"error falls only {d['amp_dn']*1e3:+.1f} mV below the "
                       f"lock point (need {MIN_ERR_AMP*1e3:.0f} mV)")
    if REQUIRE_BRACKET and not d['bracketed']:
        bad.append(f"extrema at {d['v_lo']:.4f} / {d['v_hi']:.4f} V do not "
                   f"straddle the peak at {d['v_peak']:.4f} V")
    if abs(d['slope']) < MIN_SLOPE:
        bad.append(f"|slope| {abs(d['slope']):.1f} < MIN_SLOPE {MIN_SLOPE} V/V")
    return bad


def plot_sweep(v, trans, err, d, title):
    fig, ax1 = plt.subplots(figsize=(8, 4.5))
    ax1.plot(v, trans, 'b-')
    ax1.set_xlabel('OUT2 (V)')
    ax1.set_ylabel('transmission (V)', color='b')

    ax2 = ax1.twinx()
    ax2.plot(v, err, 'r-')
    ax2.set_ylabel('error signal (V)', color='r')
    print(d)
    if d is not None:
        ax2.plot([d['v_lo'], d['v_hi']], [d['e_lo'], d['e_hi']], 'ko', ms=5)
    plt.title(title)
    plt.tight_layout()
    plt.show()




def engage(d):
    """Close the loop on the characterised resonance d."""
    if SIGN_OVERRIDE is None:
        sign = -np.sign(d['slope'])
        if sign == 0:
            raise RuntimeError("measured slope is zero - cannot set loop sign")
    else:
        sign = float(SIGN_OVERRIDE)

    set_output(d['v_peak'])
    pid.setpoint = d['e_peak'] if SETPOINT_FROM_SCAN else 0.0
    pid.p = sign * P_GAIN
    pid.i = sign * I_GAIN
    return sign


def acquire(center, span, label):
    """Sweep, rank candidates, validate, engage the best lockable one."""
    v, trans, err = sweep(center, span)
    trans_s = smooth(trans)
    
    cand = candidates(v, trans_s)
    print(f"  {len(cand)} candidate peak(s): "
          + ", ".join(f"{v[k]:.4f} V ({trans_s[k]:.3f} V)" for k in cand))

    chosen, rejected = None, []
    for k in cand:
        d = characterise(v, trans_s, err, k)
        bad = validate(d)
        tag = (f"{d['v_peak']:.4f} V  T={d['t_peak']:.3f}  "
               f"err={d['e_peak']:+.4f}  slope={d['slope']:+.1f} V/V  "
               f"disc=-{d['amp_dn']*1e3:.1f}/+{d['amp_up']*1e3:.1f} mV")
        if bad:
            print(f"  REJECT  {tag}")
            for b in bad:
                print(f"            {b}")
            rejected.append(d)
        else:
            print(f"  ACCEPT  {tag}")
            chosen = d
            break

    if chosen is None:
        print("  no resonance passed the discriminant test")
        if PLOT_ACQUISITION:
            plot_sweep(v, trans, err, rejected[0] if rejected else None,
                       label + " (ALL REJECTED)")
        return None

    if PLOT_ACQUISITION:
        plot_sweep(v, trans, err, chosen, label)

    sign = engage(chosen)
    print(f"  engaged: sign={sign:+.0f}  i={pid.i:+.3f}  p={pid.p:+.3f}  "
          f"setpoint={pid.setpoint:+.4f} V  "
          f"f_unity~{abs(pid.i)*abs(chosen['slope']):.0f} Hz")
    return chosen, (v, trans, err)
# =====================  NETWORK ANALYSER (PLANT)  ====================
def _get_na():
    """pyrpl exposes the NA under different names depending on version."""
    for owner, attr in ((p, 'networkanalyzer'), (p, 'network_analyzer'),
                        (p.rp, 'na'), (p.rp, 'networkanalyzer')):
        try:
            na = getattr(owner, attr)
            if na is not None:
                return na
        except Exception:
            pass
    return None
def _set_first(obj, names, value):
    """Set the first attribute name that exists. Returns the name used."""
    for n in names:
        if hasattr(obj, n):
            try:
                setattr(obj, n, value)
                return n
            except Exception:
                continue
    return None
def _na_frequencies(na, n):
    for attr in ('frequencies', 'data_x', 'x'):
        if hasattr(na, attr):
            f = np.asarray(getattr(na, attr))
            if f.size == n:
                return f
    # fall back to reconstructing the axis
    if NA_LOGSCALE:
        return np.logspace(np.log10(NA_START), np.log10(NA_STOP), n)
    return np.linspace(NA_START, NA_STOP, n)
def plant_sweep(label="plant"):
    """
    Measure OUT2 -> iq0 while the loop is closed at reduced gain.

    The measured response is M = P / (1 - P*C), where C = p + i/(j*f) is the
    PID. Above unity gain |P*C| << 1 so M ~ P directly; below it the servo
    suppresses the response and the de-embedded P = M / (1 + M*C) becomes
    ill-conditioned, so those points are masked out.

    The NA drive sums onto OUT2 downstream of the pid0 clamp. That is safe
    here because the RP cannot exceed +/-1 V in any case, but it does mean
    min_voltage/max_voltage do not gate the injected tone.
    """
    na = _get_na()
    if na is None:
        print("  network analyser not found on this pyrpl version - skipping")
        return None

    t0_trans, _ = measure()
    saved = dict(i=pid.i, p=pid.p, bw=iq.bandwidth)
    result = None

    try:
        print(f"\n--- {label}: measuring plant "
              f"{NA_START:.0f} Hz to {NA_STOP/1e3:.0f} kHz ---")
        # widen the demodulator so its own poles are not what we measure
        iq.bandwidth = [NA_IQ_BANDWIDTH, NA_IQ_BANDWIDTH]
        print("  iq0 bandwidth for sweep:", iq.bandwidth)
        # drop the loop gain so it is transparent over the measured band
        pid.i = np.sign(saved['i']) * NA_I_GAIN if saved['i'] != 0 else NA_I_GAIN
        pid.p = 0.0
        time.sleep(NA_SETTLE)

        t1_trans, _ = measure()
        print(f"  transmission {t0_trans:.3f} -> {t1_trans:.3f} V "
              f"at reduced gain")

        _set_first(na, ['input'], 'iq0')
        _set_first(na, ['output_direct'], 'out2')
        _set_first(na, ['start_freq', 'start'], NA_START)
        _set_first(na, ['stop_freq', 'stop'], NA_STOP)
        _set_first(na, ['points'], NA_POINTS)
        _set_first(na, ['rbw'], NA_RBW)
        _set_first(na, ['avg_per_point', 'average_per_point', 'avg'], NA_AVG)
        _set_first(na, ['amplitude'], NA_AMPLITUDE)
        _set_first(na, ['logscale'], NA_LOGSCALE)
        _set_first(na, ['acbandwidth'], 0)

        t0 = time.perf_counter()
        z = np.asarray(na.single())
        print(f"  sweep took {time.perf_counter()-t0:.1f} s, "
              f"{len(z)} points")
        f = _na_frequencies(na, len(z))

        # de-embed the servo:  P = M / (1 + M*C)
        C = saved['p'] + (np.sign(saved['i']) * NA_I_GAIN) / (1j * f)
        denom = 1.0 + z * C
        P = z / denom
        trust = np.abs(denom) > 0.2          # mask where inversion is unstable
        if not NA_DEEMBED:
            P, trust = z, np.ones_like(f, dtype=bool)

        mag = np.abs(P)
        band = trust & (f > NA_SEARCH_MIN)
        f_res, q_hint = None, None
        if band.sum() > 5:
            kk = int(np.argmax(mag[band]))
            f_res = float(f[band][kk])
            base = np.median(mag[band][:max(5, band.sum() // 10)])
            q_hint = float(mag[band][kk] / base) if base > 0 else np.nan
            # only call it a resonance if it actually stands out
            if not np.isfinite(q_hint) or q_hint < 2.0:
                f_res, q_hint = None, q_hint

        dc = float(np.median(np.abs(P)[trust & (f < 5 * NA_SEARCH_MIN)])) \
            if (trust & (f < 5 * NA_SEARCH_MIN)).sum() else np.nan

        result = dict(f=f, M=z, P=P, trust=trust, f_res=f_res,
                      peak_ratio=q_hint, dc_gain=dc, label=label)

        print(f"  |P| at low frequency: {dc:.1f} V/V "
              f"(compare with the swept discriminant slope)")
        if f_res is not None:
            f_unity = f_res / NA_TARGET_MARGIN
            print(f"  resonance at {f_res:.0f} Hz, "
                  f"{q_hint:.1f}x above the low-frequency level")
            print(f"  -> put unity gain near {f_unity:.0f} Hz")
            if dc > 0:
                print(f"  -> I_GAIN <= {f_unity/dc:.2f}")
        else:
            print("  no clear mechanical resonance in band; the limit is "
                  "elsewhere (iq filter, or slope steeper than assumed)")

        if PLOT_LOG:
            fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
            ax[0].loglog(f, np.abs(z), 'c-', lw=1, label='measured M (closed loop)')
            ax[0].loglog(f[trust], np.abs(P)[trust], 'b-', label='plant P (de-embedded)')
            ax[0].set_ylabel('|gain| (V/V)')
            ax[0].legend(fontsize=8)
            ax[0].set_title(f'{label}: OUT2 -> iq0')
            ax[1].semilogx(f[trust], np.angle(P[trust], deg=True), 'b-')
            ax[1].axhline(0, color='k', lw=0.5)
            ax[1].set_ylabel('phase (deg)')
            ax[1].set_xlabel('frequency (Hz)')
            if f_res is not None:
                for a in ax:
                    a.axvline(f_res, color='r', ls='--', lw=0.8)
                    a.axvline(f_res / NA_TARGET_MARGIN, color='g', ls=':', lw=0.8)
            plt.tight_layout()
            plt.show()

    except Exception as e:
        print("  network analyser error:", e)
    finally:
        try:
            _set_first(na, ['amplitude'], 0)
            _set_first(na, ['output_direct'], 'off')
        except Exception:
            pass
        iq.bandwidth = saved['bw']
        pid.p = saved['p']
        pid.i = saved['i']
        time.sleep(0.5)
        t2, _ = measure()
        print(f"  restored: i={pid.i:+.3f}, iq bw={iq.bandwidth}, "
              f"T={t2:.3f} V")

    return result


def shutdown():
    print("\nCleaning up Red Pitaya...")
    try:
        pid.p = 0.0
        pid.i = 0.0
        if RAMP_DOWN:
            v0 = pid.ival
            for vv in np.linspace(v0, 0.0, RAMP_DOWN_STEPS):
                pid.ival = float(vv)
                time.sleep(RAMP_DOWN_TIME / RAMP_DOWN_STEPS)
        pid.ival = 0.0
        pid.output_direct = 'off'
    except Exception as e:
        print("  pid shutdown warning:", e)
    try:
        s.stop()
    except Exception as e:
        print("  scope stop warning:", e)
    for name, mod in (("iq0", iq), ("asg1", asg)):
        try:
            mod.amplitude = 0
            mod.output_direct = 'off'
        except Exception as e:
            print(f"  {name} disable warning:", e)
    try:
        p.rp.end_all()
    except Exception as e:
        print("  close warning:", e)


 
# =============================  RUN  =================================

log_t, log_tl, log_ival, log_trans, log_raw, log_err = [], [], [], [], [], []
acq_sweep, na_start_result, na_end_result = None, None, None
lock, relocks = None, 0

try:
    print("Initial acquisition...")
    result = acquire(V_GUESS, ACQ_SPAN, "initial acquisition")
    if result is None:
        raise RuntimeError("initial acquisition failed")
    lock, acq_sweep = result
    t_peak = lock['t_peak']

    t_start = time.perf_counter()
    t_engage = t_start                       # time of the most recent engage

    if RUN_NA_AT_START:
        time.sleep(PLOT_SETTLE)              # let the creep transient pass
        na_start_result = plant_sweep("plant @ start")
        t_engage = time.perf_counter()       # gain changes count as re-engage

    print(f"\nHolding lock for {RUN_TIME} s  Terminate to stop early\n")
    cycle = 0

    while time.perf_counter() - t_start < RUN_TIME:
        time.sleep(MONITOR_DT)
        cycle += 1
        now = time.perf_counter()
        t_now = now - t_start
        trans, raw = measure()
        err = raw - pid.setpoint             # the actual PID error
        ival = pid.ival

        log_t.append(t_now)
        log_tl.append(now - t_engage)        # time since the current engage
        log_ival.append(ival)
        log_trans.append(trans)
        log_raw.append(raw)
        log_err.append(err)
        good_hist=[]

        if cycle % PRINT_EVERY == 0:
            print(f"t = {t_now:7.1f} s   out = {ival:+.4f} V   "
                  f"T = {trans:.3f} V   err = {err:+.5f} V")

        if trans > LOCK_THRESHOLD * t_peak:
            good_hist.append(ival)           # remember where resonance was
            if len(good_hist) > 10:
                good_hist.pop(0)
        else:
            relocks += 1
            home = float(np.median(good_hist)) if good_hist else V_GUESS
            print(f"  LOCK LOST (T = {trans:.3f} < "
                  f"{LOCK_THRESHOLD * t_peak:.3f}) - re-acquiring "
                  f"[{relocks}/{MAX_RELOCKS}]")
            print(f"    integrator at {ival:+.4f} V; searching around the "
                  f"last healthy point {home:.4f} V")
            if relocks > MAX_RELOCKS:
                raise RuntimeError("too many failed re-acquisitions")
            pid.p = 0.0
            pid.i = 0.0
            pid.ival = float(np.clip(home, V_MIN, V_MAX))   # undo the windup
            result = acquire(home, RELOCK_SPAN, f"relock {relocks}")
            if result is None:
                print("  narrow relock found nothing lockable - widening")
                result = acquire(home, ACQ_SPAN, f"relock {relocks} (wide)")
                if result is None:
                    raise RuntimeError("no lockable resonance found")
            lock, acq_sweep = result
            t_peak = lock['t_peak']
            good_hist = [lock['v_peak']]
            t_engage = time.perf_counter()


except KeyboardInterrupt:
    print("\nInterrupted by user.")
except Exception as e:
    print("\nERROR:", e)
finally:
    setpoint_used = pid.setpoint
    shutdown()
    timestamp = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # ---- save the accepted acquisition sweep ------------------------
    if SAVE_SWEEP and acq_sweep is not None:
        try:
            v, tr, er = acq_sweep
            fname = SAVE_DIR + r'\cavitysweep{}.txt'.format(timestamp)
            np.savetxt(fname, np.column_stack([v, tr, er]),
                       header='out_V transmission_V error_V | '
                              'F_MOD={}Hz, phase={}deg, mod_depth={}V, '
                              'iq_bw={}Hz'.format(F_MOD, PHASE, MOD_DEPTH,
                                                  IQ_BANDWIDTH))
            print("  sweep saved to", fname)
        except Exception as e:
            print("  sweep save warning:", e)

    # ---- save the network analyser data -----------------------------
    if NA_SAVE:
        for tag, res in (('start', na_start_result), ('end', na_end_result)):
            if res is None:
                continue
            try:
                fname = SAVE_DIR + r'\cavityplant_{}_{}.txt'.format(tag, timestamp)
                np.savetxt(fname,
                           np.column_stack([res['f'],
                                            res['M'].real, res['M'].imag,
                                            res['P'].real, res['P'].imag,
                                            res['trust'].astype(int)]),
                           header='f_Hz reM imM reP imP trusted | '
                                  'amp={}V, na_i_gain={}, iq_bw={}Hz, '
                                  'i_gain={}'.format(NA_AMPLITUDE, NA_I_GAIN,
                                                     NA_IQ_BANDWIDTH, I_GAIN))
                print(f"  plant ({tag}) saved to", fname)
            except Exception as e:
                print(f"  plant save warning ({tag}):", e)

    # ---- log, gated on the settled lock -----------------------------
    if len(log_t) > 0:
        log = np.column_stack([log_t, log_tl, log_ival,
                               log_trans, log_raw, log_err])
        settled = log[:, 1] >= PLOT_SETTLE

        if SAVE_LOG:
            try:
                fname = SAVE_DIR + r'\cavitylock{}.txt'.format(timestamp)
                np.savetxt(
                    fname, log,
                    header='t_s t_since_engage_s ival_V transmission_V '
                           'raw_iq0_V pid_error_V | '
                           'F_MOD={}Hz, phase={}deg, mod_depth={}V, '
                           'iq_bw={}Hz, i_gain={}, p_gain={}, setpoint={}V, '
                           'plot_settle={}s'.format(
                               F_MOD, PHASE, MOD_DEPTH, IQ_BANDWIDTH,
                               I_GAIN, P_GAIN, setpoint_used, PLOT_SETTLE))
                print("  log saved to", fname)
            except Exception as e:
                print("  log save warning:", e)

        if PLOT_LOG and settled.sum() > 2:
            # blank the unsettled samples so gaps show where it re-acquired
            show = log.copy()
            show[~settled, 2:] = np.nan
            fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 7))
            ax[0].plot(show[:, 0], show[:, 2])
            ax[0].set_ylabel('OUT2 (V)')
            ax[0].set_title(f'cavity drift, I_GAIN = {I_GAIN}  '
                            f'(first {PLOT_SETTLE:.0f} s after each engage hidden)')
            ax[1].plot(show[:, 0], show[:, 3])
            ax[1].set_ylabel('transmission (V)')
            ax[2].plot(show[:, 0], show[:, 5])
            ax[2].axhline(0, color='k', lw=0.5)
            ax[2].set_ylabel('PID error (V)')
            ax[2].set_xlabel('time (s)')
            plt.tight_layout()
            plt.show()

        # ---- report -------------------------------------------------
        print("\n" + "=" * 62)
        print("REPORT")
        print("=" * 62)
        if lock is not None:
            print(f"  transmission      {lock['t_peak']:.3f} V")
            print(f"  discriminant      {lock['slope']:+.1f} V/V, "
                  f"-{lock['amp_dn']*1e3:.1f} / +{lock['amp_up']*1e3:.1f} mV")
            print(f"  setpoint          {setpoint_used:+.4f} V")
            print(f"  I_GAIN            {I_GAIN}  "
                  f"-> f_unity ~ {I_GAIN*abs(lock['slope']):.0f} Hz")
        print(f"  re-acquisitions   {relocks}")

        if settled.sum() > 2:
            iv, tt, ee = log[settled, 2], log[settled, 3], log[settled, 5]
            ts = log[settled, 0]
            print(f"  settled samples   {settled.sum()} of {len(log)} "
                  f"({PLOT_SETTLE:.0f} s discarded after each engage)")
            print(f"  drift             {1e3*(iv.max()-iv.min()):.1f} mV p-p, "
                  f"{1e3*np.polyfit(ts, iv, 1)[0]:.3f} mV/s")
            print(f"  transmission      {tt.mean():.4f} V mean, "
                  f"{100*tt.std()/tt.mean():.2f} % rms")
            print(f"  PID error         {1e3*ee.mean():+.3f} mV mean, "
                  f"{1e3*ee.std():.3f} mV rms")
            if lock is not None and lock['slope'] != 0:
                print(f"                    = {1e6*ee.std()/abs(lock['slope']):.1f} "
                      f"uV of detuning at OUT2")
            print(f"  NOTE: sampled at {MONITOR_DT:.2f} s, so everything above "
                  f"{0.5/MONITOR_DT:.1f} Hz is aliased into these numbers.")
        else:
            print("  not enough settled samples for statistics - "
                  "lower PLOT_SETTLE or run longer")

        for tag, res in (('start', na_start_result), ('end', na_end_result)):
            if res is None:
                continue
            print(f"  plant ({tag})      |P| ~ {res['dc_gain']:.1f} V/V at low f")
            if res['f_res'] is not None:
                print(f"                    resonance {res['f_res']:.0f} Hz "
                      f"({res['peak_ratio']:.1f}x)")
                print(f"                    suggested f_unity "
                      f"{res['f_res']/NA_TARGET_MARGIN:.0f} Hz, "
                      f"I_GAIN <= {res['f_res']/NA_TARGET_MARGIN/res['dc_gain']:.2f}")
            else:
                print("                    no clear resonance in band")
        print("=" * 62)

# %%
