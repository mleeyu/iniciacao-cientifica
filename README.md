# Iniciação Científica

Este repositório contém os códigos e dados do projeto de iniciação científica "Avaliação das Ferramentas de *Backtesting* para o Valor-em-Risco e o *Expected Shortfall*".

O script principal é `main.R`, que depende de três componentes:
- O pacote `parallelDGP`;
- O script `backtests.R`;
- O script `scores.R`.

E os dados estão em `./data/output.RData`.

## Requisitos
### Dependências do `parallelDGP`
  - Biblioteca GSL (GNU Scientific Library), instalada globalmente ou localmente;
  - Pacotes R:
    - `Rcpp`;
    - `RcppEigen`;
    - `RcppGSL`;
    - `nloptr`.
### Dependências do `backtests.R`
  - Pacotes R:
    - `MASS`;
    - `quantreg`;
    - `GAS`;
    - `esback`.
