# Iniciação Científica

Este repositório contém os códigos e dados do projeto de iniciação científica "Avaliação das Ferramentas de *Backtesting* para o Valor-em-Risco e o *Expected Shortfall*".

Os scripts `monte_carlo_simulations.R` e `empirical_applications.R` contêm, respectivamente, as simulações de Monte Carlo e as aplicações empíricas.

O script `monte_carlo_simulations.R` depende do:
- Pacote `parallelDGP`;
- Script `./utils/calibration_tests.R`;
- Script `./utils/scoring_functions.R`.

O script `empirical_applications.R` depende do:
- Pacote `parallelDGP`;
- Script `./utils/calibration_tests.R`;
- Script `./utils/scoring_functions.R`;
- Script `./utils/outlier_detector.R`.

Os resultados estão disponíveis em `./data/*.RData` e as tabelas geradas com o script `./utils/table_maker.R`.

## Requisitos
### Dependências do `parallelDGP`
  - Biblioteca GSL (GNU Scientific Library), instalada globalmente ou localmente;
  - Pacotes R:
    - `Rcpp`;
    - `RcppEigen`;
    - `RcppGSL`;
    - `nloptr`.
### Dependências do `./utils/calibration_tests.R`
  - Pacotes R:
    - `MASS`;
    - `quantreg`;
    - `GAS`;
    - `esback`.
