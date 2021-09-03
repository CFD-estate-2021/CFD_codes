Ho costruito un codice per convezione-diffusione per un dominio rettangolare.
Come puoi vedere dai valori stampati nella cella Gauss-Seidel i risulatati delle iterazioni in una cella divergono.
Ho notato che le celle del contorno non sono state programmate correttamente e sto correggendo questo adesso.
Spero che così non ci siano più problemi di convergenza numerica.

Il dominio è così strutturato:
- faccia SX: T = 303 K
- faccia DX: T = 273 K
- faccia inferiore: q_out = 50 J/s
- faccia superiore: q_out = 50 J/s
