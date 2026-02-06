import re
import io
import numpy as np
import pandas as pd
import copy


import streamlit as st
import plotly.graph_objects as go

from scipy.signal import savgol_filter
from scipy.optimize import curve_fit


# -----------------------------
# Pré-processamento: outliers + smoothing
# -----------------------------
def hampel_filter_1d(x, window_size=3, n_sigmas=3.0):
    x = x.astype(float).copy()
    n = len(x)
    k = 1.4826
    for i in range(n):
        i0 = max(i - window_size, 0)
        i1 = min(i + window_size + 1, n)
        window = x[i0:i1]
        med = np.median(window)
        mad = np.median(np.abs(window - med))
        if mad == 0:
            continue
        sigma = k * mad
        if np.abs(x[i] - med) > n_sigmas * sigma:
            x[i] = med
    return x


def preprocess_spectrum(y, hampel_win=3, hampel_sig=3.0, sg_window=7, sg_poly=2):
    y1 = hampel_filter_1d(y, window_size=hampel_win, n_sigmas=hampel_sig)

    sg_window = int(sg_window)
    if sg_window % 2 == 0:
        sg_window += 1
    if sg_window >= len(y1):
        sg_window = len(y1) if (len(y1) % 2 == 1) else (len(y1) - 1)
    sg_window = max(sg_window, 5)

    sg_poly = int(sg_poly)
    sg_poly = min(sg_poly, sg_window - 2)
    sg_poly = max(sg_poly, 2)

    y2 = savgol_filter(y1, window_length=sg_window, polyorder=sg_poly, mode="interp")
    return y2


# -----------------------------
# Leitura dos arquivos
# -----------------------------
def read_calibration_xlsx(uploaded_file):
    df = pd.read_excel(uploaded_file)
    df.columns = [str(c).strip() for c in df.columns]

    wl_col = df.columns[0]
    wavelengths = df[wl_col].to_numpy(dtype=float)

    spectra_cols = df.columns[1:]
    hazen_vals = []
    spectra = []

    for c in spectra_cols:
        m = re.search(r"Hazen\s*([0-9]+)", c, re.IGNORECASE)
        if not m:
            raise ValueError(
                f"Coluna '{c}' não parece conter 'Hazen N'. "
                "Renomeie para 'Hazen 0', 'Hazen 20', etc."
            )
        haz = int(m.group(1))
        hazen_vals.append(haz)
        spectra.append(df[c].to_numpy(dtype=float))

    hazen_vals = np.array(hazen_vals, dtype=float)
    spectra = np.column_stack(spectra)

    order = np.argsort(hazen_vals)
    hazen_vals = hazen_vals[order]
    spectra = spectra[:, order]

    return wavelengths, hazen_vals, spectra


def read_samples_xlsx(uploaded_file):
    df = pd.read_excel(uploaded_file)
    df.columns = [str(c).strip() for c in df.columns]

    wl_col = df.columns[0]
    wavelengths = df[wl_col].to_numpy(dtype=float)
    sample_cols = df.columns[1:]
    Y = np.column_stack([df[c].to_numpy(dtype=float) for c in sample_cols])
    return wavelengths, sample_cols.tolist(), Y


# -----------------------------
# Índice na banda (400–430 por default)
# -----------------------------
def band_attenuation_index(spectrum, ref_spectrum):
    ref = np.clip(ref_spectrum, 1e-9, None)
    att = np.mean((ref - spectrum) / ref)
    return float(att)


# -----------------------------
# Modelos candidatos
# -----------------------------
def f_linear(x, a, b):
    return a + b * x


def f_poly2(x, a, b, c):
    return a + b * x + c * (x**2)


def f_poly3(x, a, b, c, d):
    return a + b * x + c * (x**2) + d * (x**3)


def f_log(x, a, b, c):
    return a + b * np.log(x + c)


def f_power(x, a, b, c):
    return a * np.power(x + c, b)


def loocv_rmse(x, y, func, p0, bounds):
    errs = []
    n = len(x)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        xi, yi = x[mask], y[mask]
        try:
            popt, _ = curve_fit(func, xi, yi, p0=p0, bounds=bounds, maxfev=20000)
            yhat = func(x[i], *popt)
            errs.append((yhat - y[i]) ** 2)
        except Exception:
            return np.inf
    return float(np.sqrt(np.mean(errs))) if errs else np.inf


def fit_best_curve(x, y, enabled_models):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    eps = 1e-6
    minx = float(np.min(x))

    candidates_all = {
        "linear": (f_linear, [np.mean(y), 1.0], (-np.inf, np.inf)),
        "poly2":  (f_poly2,  [np.mean(y), 1.0, 1.0], (-np.inf, np.inf)),
        "poly3":  (f_poly3,  [np.mean(y), 1.0, 1.0, 1.0], (-np.inf, np.inf)),
        "log":    (f_log,    [np.mean(y), 50.0, max(eps, -minx + eps)],
                   ([-np.inf, -np.inf, -minx + eps], [np.inf, np.inf, np.inf])),
        "power":  (f_power,  [max(np.mean(y), eps), 1.0, max(eps, -minx + eps)],
                   ([0.0, -5.0, -minx + eps], [np.inf, 5.0, np.inf])),
    }

    scores = []
    for name in enabled_models:
        func, p0, bounds = candidates_all[name]
        rmse = loocv_rmse(x, y, func, p0, bounds)
        scores.append((name, rmse, func, p0, bounds))

    scores.sort(key=lambda t: t[1])
    best_name, best_rmse, best_func, best_p0, best_bounds = scores[0]

    popt, _ = curve_fit(best_func, x, y, p0=best_p0, bounds=best_bounds, maxfev=20000)

    return {
        "best_name": best_name,
        "best_rmse_loocv": best_rmse,
        "func": best_func,
        "popt": popt,
        "all_scores": [(n, r) for (n, r, *_rest) in scores],
    }


# -----------------------------
# UI
# -----------------------------
st.set_page_config(page_title="Hazen Measurer (400–430 nm)", layout="wide")
st.title("🧪 Hazen Measurer — Banda 400–430 nm")
st.caption("Calibração automática (linear/log/power/polinomial), com remoção de outliers e suavização do espectro.")

with st.sidebar:
    st.header("⚙️ Configurações")

    st.subheader("Faixa espectral")
    wl_min = st.number_input("λ mínimo (nm)", value=400.0, step=1.0)
    wl_max = st.number_input("λ máximo (nm)", value=430.0, step=1.0)

    st.subheader("Filtro de outliers (Hampel)")
    hampel_win = st.slider("Janela (pontos p/ lado)", min_value=1, max_value=15, value=3)
    hampel_sig = st.slider("n_sigmas", min_value=1.0, max_value=6.0, value=3.0, step=0.5)

    st.subheader("Suavização (Savitzky–Golay)")
    sg_window = st.slider("Window length (ímpar)", min_value=5, max_value=51, value=7, step=2)
    sg_poly = st.slider("Poly order", min_value=2, max_value=5, value=2)

    st.subheader("Modelos para testar")
    enabled_models = st.multiselect(
        "Selecione modelos candidatos",
        options=["linear", "log", "power", "poly2", "poly3"],
        default=["linear", "log", "power", "poly2", "poly3"],
    )

    st.subheader("Referência")
    ref_strategy = st.selectbox(
        "Como escolher a referência?",
        options=["Menor Hazen (recomendado)", "Escolher manualmente"],
        index=0
    )
    ref_manual_value = None
    if ref_strategy == "Escolher manualmente":
        ref_manual_value = st.number_input("Valor Hazen de referência", value=0.0, step=5.0)

    st.divider()
    st.markdown("**Dica:** Garanta que no arquivo de calibração as colunas estejam como `Hazen 0`, `Hazen 20`, etc.")


colA, colB = st.columns([1, 1])

with colA:
    st.subheader("1) Upload — Calibração")
    calib_file = st.file_uploader("Arquivo de calibração (.xlsx)", type=["xlsx"], key="calib")

with colB:
    st.subheader("2) Upload — Amostras reais")
    sample_file = st.file_uploader("Arquivo de amostras (.xlsx)", type=["xlsx"], key="samples")

if not calib_file:
    st.info("Envie o **arquivo de calibração** para começar.")
    st.stop()

if not enabled_models:
    st.warning("Selecione pelo menos um modelo para testar.")
    st.stop()

# -----------------------------
# Processa calibração
# -----------------------------
try:
    wl, hazen, spectra = read_calibration_xlsx(calib_file)
except Exception as e:
    st.error(f"Erro ao ler calibração: {e}")
    st.stop()

band_mask = (wl >= wl_min) & (wl <= wl_max)
if band_mask.sum() < 8:
    st.error("Poucos pontos na faixa escolhida. Aumente a faixa ou confira o arquivo.")
    st.stop()

wl_b = wl[band_mask]
spectra_b = spectra[band_mask, :]

spectra_pp = np.zeros_like(spectra_b)
for j in range(spectra_b.shape[1]):
    spectra_pp[:, j] = preprocess_spectrum(
        spectra_b[:, j],
        hampel_win=hampel_win,
        hampel_sig=hampel_sig,
        sg_window=sg_window,
        sg_poly=sg_poly
    )

# referência
if ref_strategy == "Menor Hazen (recomendado)":
    ref_idx = int(np.argmin(hazen))
else:
    # acha a coluna mais próxima do valor digitado
    ref_idx = int(np.argmin(np.abs(hazen - float(ref_manual_value))))

ref = spectra_pp[:, ref_idx]
ref_hazen = float(hazen[ref_idx])

# índice por padrão
x_index = np.array(
    [band_attenuation_index(spectra_pp[:, j], ref) for j in range(spectra_pp.shape[1])],
    dtype=float
)

# ajusta curva
fit = fit_best_curve(x_index, hazen, enabled_models)

# -----------------------------
# Layout resultados calibração
# -----------------------------
st.divider()
st.subheader("📈 Curva de calibração")

left, right = st.columns([1.2, 0.8])

with left:
    x_min, x_max = float(np.min(x_index)), float(np.max(x_index))
    x_grid = np.linspace(x_min, x_max, 400)
    y_grid = fit["func"](x_grid, *fit["popt"])

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x_index, y=hazen, mode="markers",
        name="Padrões Hazen",
        hovertemplate="Índice=%{x:.6f}<br>Hazen=%{y:.2f}<extra></extra>"
    ))
    fig.add_trace(go.Scatter(
        x=x_grid, y=y_grid, mode="lines",
        name=f"Ajuste: {fit['best_name']}",
        hovertemplate="Índice=%{x:.6f}<br>Hazen=%{y:.2f}<extra></extra>"
    ))

    fig.update_layout(
        height=520,
        xaxis_title=f"Índice de atenuação relativo (λ={wl_min:.0f}–{wl_max:.0f} nm), ref=Hazen {ref_hazen:.0f}",
        yaxis_title="Hazen (APHA)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(l=10, r=10, t=10, b=10),
    )

    st.plotly_chart(fig, use_container_width=True)

with right:
    st.markdown("### ✅ Modelo escolhido")
    st.metric("Melhor modelo", fit["best_name"])
    st.metric("RMSE (LOOCV)", f"{fit['best_rmse_loocv']:.4f}")

    st.markdown("### 🏁 Ranking (LOOCV RMSE)")
    rank_df = pd.DataFrame(fit["all_scores"], columns=["Modelo", "RMSE_LOOCV"]).sort_values("RMSE_LOOCV")
    st.dataframe(rank_df, use_container_width=True, hide_index=True)

    st.markdown("### ℹ️ Referência")
    st.write(f"Referência utilizada: **Hazen {ref_hazen:.0f}**")


# -----------------------------
# Se tiver amostras, prediz e plota ponto
# -----------------------------
if sample_file:
    st.divider()
    st.subheader("🔎 Predição em amostras reais")

    try:
        wl_s, sample_names, Y = read_samples_xlsx(sample_file)
    except Exception as e:
        st.error(f"Erro ao ler amostras: {e}")
        st.stop()

    band_mask_s = (wl_s >= wl_min) & (wl_s <= wl_max)
    if band_mask_s.sum() < 8:
        st.error("Poucos pontos na faixa escolhida nas amostras. Ajuste a faixa ou confira o arquivo.")
        st.stop()

    wl_sb = wl_s[band_mask_s]
    Yb = Y[band_mask_s, :]

    # interpola para a grade da calibração se necessário
    if len(wl_sb) != len(wl_b) or np.max(np.abs(wl_sb - wl_b)) > 1e-9:
        Y_interp = np.zeros((len(wl_b), Yb.shape[1]))
        for j in range(Yb.shape[1]):
            Y_interp[:, j] = np.interp(wl_b, wl_sb, Yb[:, j])
        Yb = Y_interp

    xs, preds = [], []
    for j in range(Yb.shape[1]):
        y_pp = preprocess_spectrum(
            Yb[:, j],
            hampel_win=hampel_win,
            hampel_sig=hampel_sig,
            sg_window=sg_window,
            sg_poly=sg_poly
        )
        xi = band_attenuation_index(y_pp, ref)
        yi = float(fit["func"](xi, *fit["popt"]))
        xs.append(xi)
        preds.append(yi)

    out_df = pd.DataFrame({
        "Amostra": sample_names,
        "Indice_atenuacao": xs,
        "Hazen_predito": preds
    }).sort_values("Hazen_predito")

    st.dataframe(out_df, use_container_width=True, hide_index=True)

    # adiciona pontos das amostras no gráfico
    fig2 = copy.deepcopy(fig)

    fig2.add_trace(go.Scatter(
        x=xs, y=preds, mode="markers+text",
        text=sample_names, textposition="top center",
        name="Amostras reais (preditas)",
        marker=dict(symbol="x", size=10),
        hovertemplate="Amostra=%{text}<br>Índice=%{x:.6f}<br>Hazen=%{y:.2f}<extra></extra>"
    ))

    st.plotly_chart(fig2, use_container_width=True)

    # download CSV
    csv_bytes = out_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇️ Baixar resultados (CSV)",
        data=csv_bytes,
        file_name="hazen_predicoes.csv",
        mime="text/csv"
    )

else:
    st.info("Envie também o **arquivo de amostras reais** para o app calcular o Hazen e marcar os pontos na curva.")
