# Referencias Sanger wild-type

Dejá acá los 8 archivos `.fsa` (o `.ab1`) de las secuencias de referencia — uno por cada
combinación variante × orientación. La app los carga **automáticamente**; no hace falta
subirlos desde la interfaz en cada análisis.

## Convención de nombres (exacta)

| Variante | Forward | Reverse |
|---|---|---|
| DPYD*2A         | `DPYD_2A_ref_fwd.fsa`       | `DPYD_2A_ref_rev.fsa`       |
| DPYD*13A        | `DPYD_13_ref_fwd.fsa`       | `DPYD_13_ref_rev.fsa`       |
| c.2846A>T       | `DPYD_c2846AT_ref_fwd.fsa`  | `DPYD_c2846AT_ref_rev.fsa`  |
| HapB3 c.1236G>A | `DPYD_HapB3_ref_fwd.fsa`    | `DPYD_HapB3_ref_rev.fsa`    |

También aceptamos extensión `.ab1` con los mismos nombres.

## Cómo actualizar las referencias

Reemplazá el archivo correspondiente. La próxima corrida lo toma automáticamente.
Si querés deshabilitar una referencia temporalmente, renombrala con sufijo `.disabled`.
