# Modelos mátematicos del crecimiento del cáncer sólido 

El cáncer es un conjunto de enfermedades que se caracteriza por la presencia de células anormales que se diseminan, crecen y se dividen sin control. 
Un ser humano adulto tiene alrededor de un trillón de células en el cuerpo, algunas de ellas como las neuronas presentan una división muy lenta mientrás que otras tienen una capacidad de reproducción mucho más rápida, tal es el caso de las células epiteliales. Todo proceso de reproducción celular esta regulado por puntos de control que procuran la integridad de cada tejido en los organismos y mantienen bajo control cualquier forma anormal de crecimiento. 

Los humanos poseen cientos de diferentes clases de células con una morfología característica, una relación única con otras células y respuestas especificas ante la presencia de hormonas, productos químicos y estímulos que provienen del ambiente celular.

Las células pueden sufrir daños en su morfología y en su información genética que puede desencadenar el desarrollo de enfermedades o cambios estructurales. La modificación en la información que contiene la célula se denota principalmente por mutaciones que modifican el ciclo celular de las células normales, estas mutaciones son cambios en el ADN que pueden ser somáticas afectando solamente a las células del propio indviduo y creando lineas celulares con diferente genotipo, o bien, hereditarias que presentan estos cambios y pueden ser heredadas a la descendencia del sujeto. Cuando una célula sufre una mutación todas sus células hijas la heredan por lo que, mientras más temprano se genere la mutación en la vida del individuo mayor será la cantidad de células con esta mutación que pueden ser potencialmente malignas.


En la literatura se han reportado distintos modelos matematicos que intentan explicar la dinámica de la proliferación de células cáncerosas de diferentes tipos, ya que, el crecimiento célular de las células malignas puede variar por diferentes motivos como el tipo de mutaciones, la exposicion a tóxinas o el lugar de localización del desarrollo del cáncer, siendo este último uno de los parámetros más importantes para el análisis de la dinámica en el desarrollo del cáncer. 

Los modelos matemáticos permiten tener un mejor entendimiento sobre el el proceso dinámico de las células cancerosas, potencialmente ayudar a predecir el tamaño de los tumores y el tipo de tratamiento más efectivo según las afectaciones que tienen las distintas terapias contra el cáncer, tales como la inmunoterapia, la radioterapia y la inmunoterapia por mencionar algunos tipos de tratamientos.


Consideramos 5 diferentes tipos de modelos deterministas diseñados para predecir la tasa del crecimiento de tumores sólidos con respecto a los cambios en el tiempo  modelados mediante EDO's que estudian la evolución del volumen del tumor \textit{V} a través del tiempo t asumiendo que el tumor ya esta presente y tiene un volumen inicial V(t = 0) = V_0, además estos modelos asumen parametros positivos, es decir, que el volumen de los tumores crece y no se reduce por sí mismo.

Los cinco modelos tomados para el crecimiento son los siguientes:

- Exponencial: consta de una descripción natural de las etapas tempranas del desarrollo del cáncer donde el crecimiento de la población de células cáncerigenas es proporcional a la población del tumor, este modelo es funcional en las primera etapas antes de que se presente un agotamiento de recursos o la muerte celular por la acción del sistema inmune. Presenta una buena descripción de los carcinomas y melanomas.
                    
- Logistico: este modelo asume una disminución en la relación de crecimiento relativa al tamaño de la población y el máximo tamaño de la población esta limitado por un factor V_maxc conocido como apacidad de carga. A mostrado buenos resultados en tumores primarios del cerebro y gliomas.
                    
- Von Bertalanffy: intenta balancear la sintesis y la destrucción tomando en cuenta que el crecimiento del tumor es proporcional al área superficial del tumor mientrás llegan constantemente nutrientes y la muerte de esta celulas el proporcional al tamaño. Presento resultados positivos al modelar el crecimiento de tumores de cáncer de riñón, cáncer de páncreas y cáncer de útero 
                    
- Lineal: se asume un crecimiento constante de orden cero por la facilidad de propagación del cáncer en su lugar de desarrollo, debido a estas carcateristicas describe adecuadamente los tumores desarrollados en el carcinoma metastásico de células renales y cáncer de ojos.
                    
- Gompertz: este modelo muestra una reducción exponencial de la tasa de crecimiento. Presenta una buena aproximación en el crecimiento del cáncer de mama, cáncer de pulmón y linfomas.
                    
Según las condiciones del cáncer desarrollado se plantean diferentes estrategias de tratamiento, entre las cuales se encuentran la radioterapia, la quimioterapia y la inmunoterapia. En este programa se modelo el crecimiento de los tumores sólidos asociados a un tipo de crecimiento celular, cálculando la tasa de crecimiento tumoral a partir de datos reales del tamaño del tumor en dos momentos diferentes. 

Posteriormente, se necesitan datos de estrategias de tratamiento para evaluar como afecta un tratmiento en el crecimiento tumoral. En todos los trtamientos se contemplo la efectividad partiendo de una maxima efectividad y una mayor cantidad de muerte de células tumorales en el primer momento de su aplicación y la constante reducción de esta efectividad con el paso del tiempo, tomando estas condiciones como una forma periódica en cada ciclo de tratamiento.

Este código presenta aproximaciones teóricas que deben tomarse con cautela. El presente trbajo es una primera versión y se presentarán mejores apróximaciones, modelos y soluciones númericas en un futuro. 
