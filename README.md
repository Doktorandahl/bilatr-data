# bilatr theta trajectories: latent bilateral relations, 1991-2025

Data version **v2026.11** (released 2026-11-04), built for the R package
[bilatr](https://github.com/doktorandahl/bilatr) 1.0.0.

Yearly latent trajectories (theta) of the level of bilateral relations for
every country pair with at least one recorded event, 1991-2025, estimated
with bilatr's hierarchical dynamic item-response model from GDELT 1.0
event data. **Higher values of theta mean more cooperative relations.**

| file | contents |
|---|---|
| `theta_directed.rds`, `theta_directed.csv.gz` | directed model: 355,821 rows, 16,750 sender-receiver dyads |
| `theta_undirected.rds`, `theta_undirected.csv.gz` | undirected model: 195,106 rows, 9,139 unordered pairs |
| `manifest.json` | version, files, sizes, checksums, row and dyad counts |
| `SHA256SUMS`, `MD5SUMS` | checksums of every file |

Both tables cover the same 189 countries (see Known gaps).

## Loading the data

With bilatr (downloads once, verifies the checksum, and caches the file):

```r
library(bilatr)
theta <- bilatr_theta("directed", version = "v2026.11")
```

Without bilatr:

```r
base <- "https://github.com/doktorandahl/bilatr-data/releases/download/v2026.11/"
theta <- readRDS(url(paste0(base, "theta_directed.rds")))   # a tibble
theta <- read.csv(paste0(base, "theta_directed.csv.gz"), check.names = FALSE)
```

Outside R, any CSV reader that handles gzip will do. The interval columns
have non-syntactic names (`5%`, `95%`, ...); in R, quote them with
backticks.

## What the data are

**Events.** GDELT 1.0 events (the historical backfiles and the daily
updates) dated 1991-01-01 through 2025-12-31, accessed 2026-01-31. An event
is kept when both actors carry `GOV`, `MIL` or `SPY` in any of their
`Type1Code`/`Type2Code`/`Type3Code` fields and the two actors' country codes
differ (`bilatr::read_gdelt()` with its defaults). All events are used, root
events or not, with no filter on event codes. Events were further
restricted to Correlates of War states (see Known gaps).

**Model.** bilatr's `stable_gamma` model (a random-walk latent state per
dyad, a Dirichlet-multinomial likelihood over event classes, and
country-level category offsets), at yearly resolution, with events grouped
by bilatr's `ModifiedRootCode` scheme (16 of its 18 classes occur in the
data, so A = 16), reference category 5 ("Consult: meet, discuss, or
visit"), `min_n_events = 1`, and 4 chains of 1,000 warmup and 1,000
sampling iterations each. Because the model includes country-level
category offsets, theta is measured net of each country's baseline mix of
event types.

**Orientation and scale.** Higher theta makes the reference category,
a cooperative class, more likely relative to the other action types, so
higher values mean more cooperative relations. The sign is fixed by the
model's orientation constraint (the reference class's discrimination is
non-negative), the location by the prior on the initial states (mean
pinned at 0 across dyads), and the scale by normalizing the discrimination
parameters to root-mean-square 1.

**Coverage.** Dyad-years before a dyad's first recorded event are dropped:
the likelihood is evaluated only in years with events, so theta before the
first event carries no information from the dyad's own data. After the
first event every year through 2025 is kept, including years without
events (where theta is interpolated or, at the end of a series,
extrapolated by the random walk), so each dyad's series has no gaps and
ends in 2025. All numeric columns are rounded to 3 significant digits.

**Columns.** One row per dyad-year, sorted by `dyad` then `year`.

| column | description |
|---|---|
| `dyad` | Two three-letter country codes joined by `_`. The codes are GDELT's `Actor1CountryCode`/`Actor2CountryCode` (CAMEO country codes); every code in these tables is also the ISO 3166-1 alpha-3 code. Directed: `SENDER_RECEIVER`, so `SWE_NOR` is Sweden toward Norway. Undirected: the two codes in alphabetical (C-locale) order, so `NOR_SWE` only. |
| `year` | Integer, 1991-2025. |
| `mean` | Posterior mean of the smoothed theta (given all years). |
| `5%`, `95%` | 5th and 95th percentiles of the smoothed theta posterior (a 90% credible interval). |
| `filtered_mean` | Posterior mean of the forward-filtered theta, which uses only events up to and including that year: for each posterior draw a linear-Bayes filter is run forward from 1991, conditional on that draw's hyperparameters, and the draw-wise filtered state means are averaged. Use this when no look-ahead is allowed (for example as a predictor in forecasting). |
| `filtered_5%`, `filtered_95%` | 5th and 95th percentiles, across posterior draws, of the draw-wise filtered state mean. They reflect hyperparameter uncertainty only, not the filter's own state uncertainty, so they are much narrower than the smoothed interval and stay constant through years without events. They are not a credible interval for the filtered theta. |

## Known gaps

- **Correlates of War restriction.** Events were further restricted to
  those where both actors' country codes map to a Correlates of War state
  system member, so GDELT's regional and other non-state codes do not
  occur. This restriction was applied when the event table was built;
  `bilatr::read_gdelt()` does not apply it on its own.
- **Six missing states.** Serbia, Romania, Bosnia and Herzegovina,
  Montenegro, Slovenia and Timor-Leste (`SRB`, `ROU`, `BIH`, `MNE`, `SVN`,
  `TLS`) are absent from both tables, although they are Correlates of War
  system members. Their CAMEO country codes were not matched in the
  country-code crosswalk used to build the event table. So their dyads are
  missing, and their neighbors' other dyads include no events with them.

This will be corrected in the next data update (`v2027.01`), expected in
January 2027.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You may share and
adapt the data for any purpose, provided you give appropriate credit.

## How to cite

Please cite the data release by its DOI (see the Zenodo record for the
version you use; a DOI for all versions is in `CITATION.cff`), and the
bilatr package (`citation("bilatr")` in R). The underlying events are
GDELT's; please also cite:

- Leetaru, K. and Schrodt, P. A. (2013). "GDELT: Global Data on Events,
  Location and Tone, 1979-2012." Paper presented at the International
  Studies Association annual convention, San Francisco.
- The GDELT Project. Global Database of Events, Language, and Tone (GDELT).
  Accessed 31 January 2026, <https://gdeltproject.org/>.
