# Admin-Proyectos-Software-BFOA
Repositorio del proyecto _"Mejora de algoritmo BFOA"_ para la asignación _*Administración de Proyectos de Software*_, archivo main --> BFOA.py. 
<br>
<br>
## Resumen
Este proyecto consiguio optimizar el algoritmo BFOA al restructurizar el codigo y agregando la clase quimiotaxis que ayudo a conseguir mejores resultados (mejor Fitness), además de eso se agregaron algunas mejoras dentro del codigo que permitieron que el codigo se ejecutara en mucho menor tiempo y con menor numero de secuencias analizadas (NFE).
<br>
<br>
## Cambios y mejoras realizadas.
Los cambios y mejoras al codigo realizados para el algoritmo fueron las siguientes:
<br>
#### 1.- Restructuración del proyecto y inclusión de la clase "Quimiotaxis".
El cambio mas drastica que sufrio el algoritmo fue la restructuración del codigo y la inclusión de la clase quimiotaxis; la clase quimiotaxis aporto al algoritmo el permitir que las bacterias interactuen entre sí por medio de atracción y repulsión, mejorando la exploración de posibilidades y evitando soluciones poco optimas. La restructuración del codigo vino junto a la inclusión de esta clase al buscar una forma de tener el codigo mas entendible.
<br>
<br>
#### 2.- Adaptación dinámica de parámetros.
Se incluyo en el algoritmo una función lineal de adaptación dinámica de parámetros (especificamente de los de atracción y repulsión), que permite que mientras más vaya avanzando el algoritmo en sus iteraciones menos sea la cantidad de soluciones exploradas.
<br>
<br>
#### 3.- Mecanismo de mutación adaptativa.
Ahora la inclusión de gaps se realiza de manera mas inteligente al hacer que las bacterías con menor fitness tengan mayor tasa de mutación.
<br>
<br>
#### 4.- Uso de información historica.
Se almacena la mejor solución actual (mejor fitness) en una variable la cuál ayudara a guiar el proceso de Tumbo al introducir una pequeña probabilidad de que la mutación se vea influenciada por la posición de gaps encontradas en la mejor solución actual.