# Modelos mátematicos del crecimiento del cáncer sólido 

El cáncer es un conjunto de enfermedades que se caracteriza por la presencia de células anormales que se diseminan, crecen y se dividen sin control. 
Un ser humano adulto tiene alrededor de un trillón de células en el cuerpo, algunas de ellas como las neuronas presentan una división muy lenta mientrás que otras tienen una capacidad de reproducción mucho más rápida, tal es el caso de las células epiteliales. Todo proceso de reproducción celular esta regulado por puntos de control que procuran la integridad de cada tejido en los organismos y mantienen bajo control cualquier forma anormal de crecimiento. 

Los humanos poseen cientos de diferentes clases de células con una morfología característica, una relación única con otras células y respuestas especificas ante la presencia de hormonas, productos químicos y estímulos que provienen del ambiente celular.

Las células pueden sufrir daños en su morfología y en su información genética que puede desencadenar el desarrollo de enfermedades o cambios estructurales. La modificación en la información que contiene la célula se denota principalmente por mutaciones que modifican el ciclo celular de las células normales, estas mutaciones son cambios en el ADN que pueden ser somáticas afectando solamente a las células del propio indviduo y creando lineas celulares con diferente genotipo, o bien, hereditarias que presentan estos cambios y pueden ser heredadas a la descendencia del sujeto. Cuando una célula sufre una mutación todas sus células hijas la heredan por lo que, mientras más temprano se genere la mutación en la vida del individuo mayor será la cantidad de células con esta mutación que pueden ser potencialmente malignas.


En la literatura se han reportado distintos modelos matematicos que intentan explicar la dinámica de la proliferación de células cáncerosas de diferentes tipos, ya que, el crecimiento célular de las células malignas puede variar por diferentes motivos como el tipo de mutaciones, la exposicion a tóxinas o el lugar de localización del desarrollo del cáncer, siendo este uno de los parámetros más importantes para el análisis de la dinámica en el desarrollo del cáncer. 

Comunmente el cáncer inicia por la generación de mutaciones geneticas  que provocan un ritmo de crecimiento célular acelerado y sin control, los modelos matemáticos permiten tener un mejor entendimiento sobre el el proceso dinámico de las células cancerosas, ayudar a predecir el tamaño de los tumores y el tipo de tratamiento más efectivo según las afectaciones que tienen las distintas terapias contra el cáncer, tales como la inmunoterapia, la radioterapia y la inmunoterapia por mencionar algunos \cite{2}.


Existen 5 diferentes tipos de modelos deterministas diseñados para predecir la tasa del crecimiento de tumores sólidos con respecto a los cambios en el tiempo \cite{8} \cite{5} \cite{7} \cite{9} modelados mediante EDO's que estudian la evolución del volumen del tumor \textit{V} a través del tiempo \textit{t} asumiendo que el tumor ya esta presente y tiene un volumen inicial \textit{V(t)} = \textit{$V_{0}$}, además estos modelos asumen parametros positivos, es decir, que el volumen de los tumores crece y no se reduce por sí mismo.

\begin{itemize}
\item Exponencial: consta de una descripción natural de las etapas tempranas del desarrollo del cáncer donde el crecimiento de la población de células cáncerigenas es proporcional a la población del tumor, este modelo es funcional en las primera etapas antes de que se presente un agotamiento de recursos o la muerte celular por la acción del sistema inmune
\begin{equation} \label{eqn:exponencial}
    \frac{dV}{dt} = k_g V
\end{equation}

 \begin{table}[ht]
            \begin{center}
                \begin{tabular}{| | c | |}
                \hline
                \hline
                Tipo de crecimiento característico de \cite{e10}:\\ %\cite{exponencial}:
                
                Carcinomas\\ Melanomas\\ 
                    
                \hline
                \hline
                \end{tabular}
                \label{tab:fruta}
            \end{center}
        \end{table}




\item Logistico: este modelo asume una disminución en la relación de crecimiento relativa al tamaño de la población y el máximo tamaño de la población esta limitado por un factor \textit{$V_{max}$} denominado capacidad de carga
\begin{equation} \label{eqn:logistico}
    \frac{dV}{dt} = k_g V (1 - \frac{V}{V_{max}})
\end{equation}

 \begin{table}[ht]
            \begin{center}
                \begin{tabular}{| | c | |}
                \hline
                \hline
                Tipo de crecimiento característico de \cite{e11}:\\ Tumores primarios del cerebro\\ Gliomas\\ 
                    
                \hline
                \hline
                \end{tabular}
                \label{tab:fruta}
            \end{center}
        \end{table}


\item Von Bertalanffy: intenta balancear la sintesis y la destrucción tomando que el crecimiento del tumor es proporcional al area superficial del tumor mientrás llegan constantemente nutrientes y la muerte de esta celulas el proporcional al tamaño \ref{egn: von}

\begin{equation}\label{egn: von}
    \frac{dV}{dt} = k_g V ^{2/3}- b V
\end{equation}

\begin{table}[ht]
            \begin{center}
                \begin{tabular}{| | c | |}
                \hline
                \hline
                Tipo de crecimiento caracteristico de \cite{m12}:\\ Cáncer de riñón\\ Cáncer de páncreas \\ Cáncer de útero\\ 
                    
                \hline
                \hline
                \end{tabular}
                \label{tab:frut}
            \end{center}
        \end{table}

\clearpage
\item Lineal: 
Se asume un crecimiento constante de orden cero por la facilidad de propagación del cáncer en su lugar de propagación.
\begin{equation} \label{eqn: lineal}
    \frac{dT}{dt} = k_g
\end{equation}

  \begin{table}[ht]
            \begin{center}
                \begin{tabular}{| | c | |}
                \hline
                \hline
                Tipo de crecimiento característico de:\\ Carcinoma metastásico de células renales \cite{8}\\ Cáncer de ojos (primeras etapas)\cite{e9}\\ 
                    
                \hline
                \hline
                \end{tabular}
                \label{tab:fra}
            \end{center}
        \end{table}

\item Gompertz: este modelo muestra una reducción exponencial de la tasa de crecimiento
\begin{equation} \label{eqn: gom1}
    \frac{dV}{dt} = r(t)  V(t)
\end{equation}
\begin{equation} \label{eqn: gom2}
    \frac{dr}{dt} = -b r(t)
\end{equation}

\begin{table}[ht]
            \begin{center}
                \begin{tabular}{| | c | |}
                \hline
                \hline
                Tipo de crecimiento caracteristico de \cite{8}:\\ Cáncer de mama\\ Cáncer de pulmón\\ Linfomas\\
                    
                \hline
                \hline
                \end{tabular}
                \label{tab:fruta}
            \end{center}
        \end{table}

\end{itemize}


  \subsection{Radioterapia}

Según las condiciones del cáncer desarrollado se plantean estrategias de radioterapia (\textit{D(t)}) con un periódo $\omega$ y un tiempo de exposición a la radiación (\textit{L}) que podemos identificar con un modelo matemático como el siguiente \cite{e13}:
        \begin{equation}
                \frac{dV}{dt} = F(t)- a V - b D V 
        \end{equation}
        \begin{equation}
                \frac{dR}{dt} = - c R + e ND + w(t) 
        \end{equation}
donde $R$ representa la dinámica de la radiación administrada en Gray y $V$ la dinámica del crecimiento de células tumorales y F(t) es el tipo de crecimiento celular. Además, la función w(t) no sindica la dinámica del tratamiento como una función a trozos\\
\begin{center}

