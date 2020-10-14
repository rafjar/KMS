### Uwagi do programu 1:
* Najszybsza jest wersja program1-numba (wymaga modułu *numba*), natomiast daje dziwne wyniki - nie ufam jej
* Druga najszybsza wersja, to program1-threads. Natomiast jest to prymitywny podział na wątki i nie testowałem poprawności
* Najpewniejszą wersją jest program1, lecz jest najwolniejszy (dla n=4 wykonuje się ok. 1:20min, dla n=5 ok. 6min)