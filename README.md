# 🧪 Hazen Measurer

![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![Streamlit](https://img.shields.io/badge/app-Streamlit-FF4B4B)
![License](https://img.shields.io/badge/license-MIT-green)

> **EN:** A Streamlit app that estimates the **Hazen (APHA / Pt-Co) color** of liquid samples
> from UV-Vis spectra in the 400–430 nm band. It cleans the spectra (Hampel outlier filter +
> Savitzky–Golay smoothing), builds an attenuation index, picks the best calibration curve among
> linear / log / power / polynomial models by **leave-one-out cross-validation**, and predicts
> new samples, with interactive Plotly charts and CSV export.

App web para **calibrar e medir a cor Hazen (APHA)** de amostras líquidas a partir de espectros
de UV-Vis, sem planilhas manuais: você envia os padrões, o app escolhe o melhor modelo e mede as
amostras reais.

## Como funciona

```
espectros .xlsx ─► recorte da banda (400–430 nm)
                ─► filtro de Hampel (remove spikes) ─► Savitzky–Golay (suaviza)
                ─► índice de atenuação relativo ao padrão de menor Hazen
                ─► ajuste de 5 modelos ─► ranking por RMSE em LOOCV ─► melhor curva
                ─► predição das amostras + gráfico + CSV
```

- **Pré-processamento robusto:** o filtro de Hampel troca pontos fora de *n*·σ (estimado pelo MAD)
  pela mediana local; depois o Savitzky–Golay suaviza sem achatar a banda.
- **Índice de atenuação:** média de `(ref − amostra) / ref` na banda, onde `ref` é o padrão de menor Hazen.
- **Seleção de modelo:** linear, logarítmico, potência, polinomial de 2º e 3º grau, ajustados com
  `scipy.optimize.curve_fit` e comparados por *leave-one-out cross-validation*, a escolha certa
  quando há poucos padrões de calibração.
- **Interface:** todos os parâmetros (faixa espectral, janelas dos filtros, modelos candidatos,
  referência) ajustáveis na barra lateral.

## Resultado com os dados de exemplo

Com os arquivos de [`examples/`](examples) (espectros **sintéticos**, com ruído e spikes de propósito):

| Modelo | RMSE (LOOCV, Hazen) |
|---|---|
| **poly2** (escolhido) | **0,65** |
| poly3 | 1,98 |
| power | 2,87 |
| linear | 6,98 |
| log | 7,01 |

| Amostra | Hazen real | Hazen predito |
|---|---|---|
| Amostra A | 15 | 14,2 |
| Amostra B | 55 | 54,1 |
| Amostra C | 120 | 119,4 |
| Amostra D | 240 | 239,7 |

## Rodando

```bash
pip install -r requirements.txt
streamlit run app.py
```

Ou com Docker:

```bash
docker build -t hazen-measurer .
docker run -p 8501:8501 hazen-measurer
```

Abra <http://localhost:8501> e envie `examples/calibracao_exemplo.xlsx` e `examples/amostras_exemplo.xlsx`.

## Formato dos arquivos

**Calibração:** primeira coluna com o comprimento de onda (nm), uma coluna por padrão com o nome `Hazen N`:

| Wavelength (nm) | Hazen 0 | Hazen 20 | Hazen 100 | … |
|---|---|---|---|---|
| 400 | 99,8 | 97,1 | 86,4 | … |

**Amostras:** primeira coluna com o comprimento de onda, uma coluna por amostra (qualquer nome).
Se a grade de comprimentos de onda for diferente da calibração, o app interpola sozinho.

## Stack

Python · Streamlit · NumPy · SciPy · pandas · Plotly · Docker

## Licença

MIT
