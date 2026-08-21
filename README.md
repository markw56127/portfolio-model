# Portfolio Management Project

Portfolio construction and evaluation on 30 stocks plus the S&P 500, using
monthly adjusted close prices from Yahoo Finance.

## Course roadmap

This repository covers the portfolio-management half of the course. The full
sequence of topics:

1. Combining individual stocks into portfolios (risk and expected return of a portfolio)
2. Maximizing return given risk, or minimizing risk given return
3. Properties of the minimum variance set (efficient frontier)
4. The single index model (with and without short sales allowed)
5. Constant correlation model (with and without short sales allowed)
6. Multigroup and multi-index models (short sales allowed)
7. Portfolio performance

Work so far: **Project 1** covers topic 1 and the entry point to topic 2.
**Project 2** covers topics 2 and 3 — the analytic efficient frontier from
Merton (1972), "An Analytic Derivation of the Efficient Portfolio Frontier"
(JFQA 7(4)), and the mutual fund theorem. **Project 3** stays on topics 2 and 3
but adds a risk free asset: the tangency portfolio, the capital allocation
line, and tracing the frontier from two tangency portfolios.

Project 3 works on a different, smaller dataset than Projects 1 and 2 — see
Data below.

## Stock selection

Five sectors, six stocks each, giving 30 stocks. Sector labels follow Yahoo
Finance's taxonomy. The columns in the data files appear in this order, so each
sector is a contiguous block of six — convenient for the multigroup models in
topic 6.

| Sector | Tickers |
|---|---|
| Technology | NVDA, AAPL, MSFT, AVGO, AMD, NOW |
| Communication Services | GOOG, META, NFLX, VZ, DIS, T |
| Financial Services | V, JPM, MS, GS, BAC, SCHW |
| Consumer Cyclical | AMZN, TSLA, MCD, SBUX, HD, EBAY |
| Healthcare | LLY, NVO, UNH, JNJ, AMGN, PFE |

`^GSPC` (the S&P 500) is carried alongside as the market index, for a total of
31 assets.

## Data

Monthly adjusted close prices downloaded through the Shiny app in `app.R`.

| File | Period | Rows |
|---|---|---|
| `stockData_train.csv` | 2017-01 – 2021-12 | 60 months |
| `stockData_test.csv` | 2022-01 – 2026-07 | 55 months |

Projects 1 and 2 use the training file only — the first five years, as the
assignment requires. The test period is held back for out-of-sample evaluation
in topic 7.

Each file is laid out as a row index, a `Date` column, then the 31 assets, so
the assignment's starter snippet indexes correctly:

```r
a <- read.csv("stockData_train.csv", sep=",", header=TRUE)
r <- (a[-1,3:ncol(a)]-a[-nrow(a),3:ncol(a)])/a[-nrow(a),3:ncol(a)]
```

## Progress

Everything below is produced by `Project1/portfolio.R`. The rendered write-up
is `Project1/Submission/report.pdf`, knitted from `Project1/report.Rmd`.

**(b) Prices to returns.** Simple (arithmetic) returns,
`R_t = P_t / P_{t-1} - 1`. Simple rather than log returns because portfolio
return is a linear combination of asset returns only under this definition.
60 monthly prices yield 59 monthly returns across 31 assets.

**(c) Means, standard deviations, covariance matrix.** Per-asset statistics for
all 31 assets in `Project1/asset_stats.csv`; the 31×31 variance-covariance
matrix in `Project1/covmat.csv`. Portfolio construction below uses the 30×30
submatrix that excludes `^GSPC`.

**(d) Risk–return plot.** All 31 assets in mean/standard-deviation space, with
the two portfolios overlaid — `Project1/risk_return.png`.

**(e) Equal allocation portfolio.** Equal weights across the 30 stocks
(w = 1/30 each).

**(f) Minimum risk portfolio.** Weights from the closed form

```
x = (Σ⁻¹ 1) / (1' Σ⁻¹ 1)
```

where Σ is the 30×30 covariance matrix of stock returns and 1 is a 30×1 vector
of ones.

### Results (monthly, 2017-02 – 2021-12)

| | Mean | Std. dev. |
|---|---|---|
| S&P 500 | 1.36% | 4.48% |
| Equal allocation (1/30) | 2.34% | 4.73% |
| Minimum risk | 1.59% | 2.07% |

![Risk and return](Project1/risk_return.png)

Two observations worth carrying forward:

- **Diversification is visible in the numbers.** The average individual stock
  has a standard deviation near 8%, but equally weighting those same 30 stocks
  gives 4.73%. The minimum-risk portfolio reaches 2.07% — below every single
  asset in the set (the lowest individual is VZ at 4.36%).
- **The minimum-risk weights are unstable.** The solution places 0.84 on JPM
  and takes 14 short positions. This is expected behaviour for unconstrained
  minimum-variance optimization, amplified here because 59 observations are
  estimating a 30-asset covariance matrix — a T/N ratio of 1.97, where you want
  T ≫ N. The matrix does invert (condition number 831), so the math is sound,
  but the weights are partly fit to estimation noise and should not be expected
  to hold out of sample. Topics 4–6 introduce the structured covariance models
  (single index, constant correlation, multigroup) that address exactly this.

## Layout

Shared inputs live at the repo root because later projects reuse them; each
project's own code and output lives in its own folder.

| Path | Contents |
|---|---|
| `app.R` | Shiny app used to download prices from Yahoo Finance |
| `stockData_train.csv` / `stockData_test.csv` | Adjusted close prices, shared across projects |
| `renv.lock`, `.Rprofile`, `renv/` | Pinned project library |
| `Project1/report.Rmd` | Write-up, knitted to pdf |
| `Project1/portfolio.R` | Analysis for parts (b)–(f); sourced by the report |
| `Project1/returns_train.csv` | 59 × 31 monthly returns |
| `Project1/asset_stats.csv` | Per-asset mean and standard deviation |
| `Project1/covmat.csv` | 31 × 31 variance-covariance matrix |
| `Project1/risk_return.png` | Risk–return plot |
| `Project1/Submission/` | The rendered deliverable |
| `Project2/report.Rmd` | Write-up, knitted to pdf |
| `Project2/frontier.R` | Analysis for parts (a)–(h); sourced by the report |
| `Project2/fig_*.png` | Parabola, hyperbola, and mutual-fund plots |
| `Project2/Submission/` | The rendered deliverable |
| `Project3/statc183c283_5stocks.txt` | Close monthly prices, 5 stocks, 1986-01 – 2003-12 |
| `Project3/statc183c283_abc.txt` | 2,499 weight combinations summing to 1, for the cloud in (e) |
| `Project3/report.Rmd` | Write-up, knitted to pdf |
| `Project3/tangency.R` | Analysis for parts (a)–(i); sourced by the report |
| `Project3/fig_*.png` | Two-stock possibilities curve, CAL, and traced frontier |
| `Project3/Submission/` | The rendered deliverable |

Each project recomputes the returns from the shared csvs rather than importing
the other's objects, so the two stay independent. They are verifiably the same
numbers: the price matrix, the 59 × 30 return matrix, the mean vector and the
covariance matrix are bitwise identical between the projects, and the minimum
risk portfolio agrees to ten decimals though Project 1 reaches it via
$\Sigma^{-1}\mathbf{1}/\mathbf{1}'\Sigma^{-1}\mathbf{1}$ and Project 2 via
$A/C$ and $\sqrt{1/C}$.

## Reproducing

Dependencies are pinned with [renv](https://rstudio.github.io/renv/), so
packages install into a project-local library rather than the global one.

```r
renv::restore()     # installs the exact versions in renv.lock
```

Then, from the repo root:

```
Rscript Project1/portfolio.R
```

Or open `Project1/report.Rmd` in RStudio and knit, which runs the same analysis
and renders the write-up. Starting R in the repo root activates the project
library automatically via `.Rprofile` — there is no environment to activate by
hand. Note that the analysis scripts `setwd()` to their own folder, but
`rmarkdown::render()` must be called from the repo root, otherwise `.Rprofile`
never runs and the project library is not on the search path.

Knitting to pdf needs pandoc and a LaTeX installation on top of the R packages.
RStudio bundles pandoc; from a plain `Rscript` call, point at it with
`RSTUDIO_PANDOC`. BasicTeX does not ship every package the rmarkdown pdf
template loads — `framed`, `xurl` and `titling` were added to the user tree
with `tlmgr init-usertree && tlmgr --usermode install framed xurl titling`,
which needs no sudo.
