# bilatr theta trajectories: latent bilateral relations, 1991-2025

[![DOI](https://zenodo.org/badge/1445654273.svg)](https://doi.org/10.5281/zenodo.23951818)

Data version **v2026.11.1** (released 2026-11-10), built for the R package
[bilatr](https://github.com/doktorandahl/bilatr) 1.0.0.

Yearly latent trajectories (theta) of the level of bilateral relations for
every country pair with at least one recorded event, 1991-2025, estimated
with bilatr's hierarchical dynamic item-response model from GDELT 1.0
event data. **Higher values of theta mean more cooperative relations.**

| file | contents |
|---|---|
| `theta_directed.rds`, `theta_directed.csv.gz` | directed model: 365,014 rows, 17,125 sender-receiver dyads |
| `theta_undirected.rds`, `theta_undirected.csv.gz` | undirected model: 200,061 rows, 9,342 unordered pairs |
| `manifest.json` | version, files, sizes, checksums, row and dyad counts |
| `SHA256SUMS`, `MD5SUMS` | checksums of every file |

Both tables cover the same 192 countries (see Known gaps).

## Loading the data

With bilatr (downloads once, verifies the checksum, and caches the file):

```r
library(bilatr)
theta <- download_bilatr_theta("directed", version = "v2026.11.1")
```

Without bilatr:

```r
base <- "https://github.com/doktorandahl/bilatr-data/releases/download/v2026.11.1/"
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
differ (`bilatr::read_gdelt()` with its defaults; an empty country code is
filled from the actor code's prefix with `recover_country_codes = TRUE`).
All events are used, root events or not, except those with a Goldstein
score between -1 and 1 (CAMEO codes 010-015, 017, 024x, 025x, 040 and
041) and those with non-CAMEO event codes. Events are further restricted
to those in which both countries are members of the Correlates of War
state system (v2024) on the date of the event
(`bilatr::filter_gdelt_states()`), which removes GDELT's regional and other
non-state codes.

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

**Coverage.** The models are fitted on each dyad's full 1991-2025 series,
with the likelihood evaluated only in dyad-years with at least one event.
Dyad-years before a dyad's first recorded event are then dropped, since
theta before the first event carries no information from the dyad's own
data, and so are dyad-years in which either country was not a member of
the state system (for example, dyads with Czechoslovakia end in 1992, and
those with Czechia or Slovakia start in 1993 at the earliest). From the
first event on, every year is kept, including years without events (where
theta is interpolated or, at the end of a series, extrapolated by the
random walk), so each dyad's series has no gaps.
The tables have 365,014 directed and 200,061 undirected dyad-years, of which 130,848 and 76,821 have at least one observed event.
All numeric columns are rounded to 4 decimal places.

**Columns.** One row per dyad-year, sorted by `dyad` then `year`.

| column | description |
|---|---|
| `dyad` | Two three-letter country codes joined by `_`. GDELT's `Actor1CountryCode`/`Actor2CountryCode` (CAMEO country codes), harmonized to ISO 3166 alpha-3 (`bilatr::gdelt_countries`): `ROU` for Romania, `SRB` for Serbia (Correlates of War state 345, including the Federal Republic of Yugoslavia and Serbia and Montenegro, in every year), `TLS` for Timor-Leste, and `CSK` (the former ISO code) for Czechoslovakia in 1991-1992. Directed: `SENDER_RECEIVER`, so `SWE_NOR` is Sweden toward Norway. Undirected: the two codes in alphabetical (C-locale) order, so `NOR_SWE` only. |
| `year` | Integer, 1991-2025. |
| `mean`, `median` | Posterior mean and median of the smoothed theta (given all years). |
| `5%`, `95%` | 5th and 95th percentiles of the smoothed theta posterior (a 90% credible interval). |
| `filtered_mean`, `filtered_median` | Mean and median, across posterior draws, of the forward-filtered theta, which uses only events up to and including that year: for each posterior draw a linear-Bayes filter is run forward from 1991, conditional on that draw's hyperparameters, and the draw-wise filtered state means are summarized. Use these when no look-ahead is allowed (for example as a predictor in forecasting). |
| `filtered_5%`, `filtered_95%` | 5th and 95th percentiles, across posterior draws, of the draw-wise filtered state mean. They reflect hyperparameter uncertainty only, not the filter's own state uncertainty, so they are much narrower than the smoothed interval and stay constant through years without events. They are not a credible interval for the filtered theta. |

## Known gaps

- **Four absent states.** GDELT 1.0 contains no events in which the
  government actors of Slovenia, Bosnia and Herzegovina, Kosovo or
  Montenegro are the source or target, once the actor filter is applied:
  GDELT 1.0 has no country code for them. These four Correlates of War
  states are therefore absent from both tables. Their absence reflects the
  coverage of the data source, not an absence of interaction, and the
  tables say nothing about their bilateral relations.
- **Missing and incomplete days.** For 24 days, most of them between
  14 June and 1 July 2025, GDELT provides no daily event files. For a
  further 13 days in 1991-1995, the files contain events, but none that pass
  the filters. Because events are aggregated to dyad-years, these gaps
  reduce the number of events in the affected years but do not create
  missing dyad-years. The model conditions on the total number of events and
  measures relations from their composition, so a gap reduces the precision
  of the estimates for the affected years rather than shifting their level,
  provided the missing days are not systematically different from the rest
  of the year. The longest gap, 14 June to 1 July 2025, is the exception:
  it coincides with the war between Israel and Iran in June 2025, so the
  2025 estimates for the dyads involved in that conflict are based on a
  composition of events from which the war itself is largely missing.

## Changes from v2026.11

This release replaces `v2026.11`, which should no longer be used.

- **Country coverage.** `v2026.11` lacked six Correlates of War states
  (Serbia, Romania, Bosnia and Herzegovina, Montenegro, Slovenia and
  Timor-Leste) because of an error in the country-code crosswalk. The
  event table is now rebuilt from raw GDELT with `bilatr::read_gdelt()` and
  `bilatr::filter_gdelt_states()`, which recovers Serbia, Romania and
  Timor-Leste. Only the four states listed under Known gaps remain absent.
  The tables now cover 192 countries.
- **Country codes** are harmonized to ISO 3166 alpha-3, and a country is
  included only in the years it is a member of the Correlates of War state
  system (see `dyad` under Columns).
- **Event codes.** Events with a Goldstein score between -1 and 1 and events
  with non-CAMEO codes were excluded in `v2026.11` too, although its README
  said there was no filter on event codes; the filter is now part of the
  scripted pipeline and documented (see Events).
- **Columns.** New `median` and `filtered_median` columns; all numeric
  columns are rounded to 4 decimal places (previously 3 significant
  digits).

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You may share and
adapt the data for any purpose, provided you give appropriate credit.

## How to cite

Please cite the data release by the DOI of the version you use:

- `v2026.11.1`: see the Zenodo record (the DOI is minted on release)
- all versions: [doi:10.5281/zenodo.23951818](https://doi.org/10.5281/zenodo.23951818)

Randahl, D., Lundgren, M. and Vikberg, C. (2026). bilatr theta
trajectories: latent bilateral relations, 1991-2025. Data release
v2026.11.1. Zenodo. <https://doi.org/10.5281/zenodo.23951818>

Please also cite the bilatr package (`citation("bilatr")` in R). The
underlying events are GDELT's; please also cite:

- Leetaru, K. and Schrodt, P. A. (2013). "GDELT: Global Data on Events,
  Location and Tone, 1979-2012." Paper presented at the International
  Studies Association annual convention, San Francisco.
- The GDELT Project. Global Database of Events, Language, and Tone (GDELT).
  Accessed 31 January 2026, <https://gdeltproject.org/>.
