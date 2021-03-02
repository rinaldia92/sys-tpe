### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 09062294-3e5f-11eb-176f-dfcbf841f111
begin
	import Pkg
	Pkg.add.(["DSP", "FFTW", "FileIO", "WAV", "PlutoUI", "Plotly"])
	Pkg.add.(["LibSndFile", "SampledSignals", "Statistics"])
	Pkg.add(["ImageMagick", "Images", "ImageIO", "Plots", "IterTools"])
	using Plots
		plotly()
	using DSP, FFTW
	using DSP.Windows: hamming, rect, hanning
	using PlutoUI
	using WAV
	using FileIO
	using Images
	using LibSndFile
	using Statistics
	using SampledSignals
	
	
	"Pad vector with zeros on the right until its length is `n`"
	padright(x, n) = copyto!(zeros(eltype(x), n), x)
	
	"Plot discrete functions"
	stem(args...; kwargs...) = sticks(args...;
									marker=:circle,
									leg=false,
									kwargs...)

	stem!(args...; kwargs...) = sticks!(args...;
									marker=:circle,
									leg=false,
									kwargs...)

	"""
	Función módulo pero con offset (opcional)
	Manda a `t` al intervalo [from, from+length)
	sumándole o restándole algún múltiplo de `len`
	"""
	cshift(t, len, from=0) = mod(t - from, len) + from
	

	# Espectrograma
	using IterTools
	function stft(x; overlap, window, nfft, rest...)
	  nwin = length(window)
	  @assert overlap < nwin

	  res = [ fft(padright(xseg .* window, nfft))
		for xseg in partition(x, nwin, nwin - overlap)]

	  return [ res[i][j] for j in 1:nfft, i in eachindex(res)]
	end


	specplot(x::AbstractMatrix; kwargs...) = 
		@error "You are entering a Matrix (2D Array). I need a Vector (1D Array)."
	function specplot(x::AbstractVector;
		  fs=1,
		  onesided=false,
		  xaxis="Tiempo (s)",
		  yaxis="Frecuencia (Hz)",
		  window=hamming(div(length(x), 16)),
		  overlap=0.5,
		  nfft=length(window),
		  kws...)

		window isa Integer && (window = rect(window))
		overlap isa AbstractFloat && (overlap = round(Int, length(window) * overlap))

		mat = stft(x; overlap=overlap, window=window, nfft=nfft)

		fmax = fs
		if onesided
		  mat = mat[1:div(size(mat, 1) + 2, 2), :]
		  fmax = fs/2
		end

	  toffset = length(window) / 2fs
	  times = range(toffset; length=size(mat, 2), stop=length(x)/fs - toffset)
	  freqs = range(0; length=size(mat, 1), stop=fmax)

		# Reubico las frecuencias negativas arriba de todo
	  if !onesided
		freqs = cshift.(freqs, fs, -fs/2)
		ord   = sortperm(freqs)
		mat   = mat[ord, :]
		freqs = freqs[ord]
	  end

		return heatmap(times, freqs, log.(abs.(mat) .+ eps());
			  xaxis=xaxis, yaxis=yaxis,
			  seriescolor=:bluesreds, legend=true, kws...)
	 return times, freqs, mat
	end

	function specplot(x :: AbstractVector{<:AbstractFloat}; kws...)
		return specplot(convert.(Complex, x); onesided=true, kws...)
	end
	
	# Polos y ceros
	
	zeropolegain(pr) = DSP.ZeroPoleGain(pr)
	zeropolegain(z, p, g) = DSP.ZeroPoleGain(z, p, g)
	polynomialratio(zpg) = DSP.PolynomialRatio(zpg)
	function polynomialratio(b, a)
	  n = max(length(a), length(b))
	  return DSP.PolynomialRatio(padright(b, n), padright(a, n))
	end
	getpoles(zpg) = DSP.ZeroPoleGain(zpg).p
	getzeros(zpg) = DSP.ZeroPoleGain(zpg).z
	getgain(zpg) = DSP.ZeroPoleGain(zpg).k
	getnumcoefs(pr) = trimlastzeros!(reverse(DSP.PolynomialRatio(pr).b.coeffs))
	getdencoefs(pr) = trimlastzeros!(reverse(DSP.PolynomialRatio(pr).a.coeffs))
	function trimlastzeros!(a)
	  !iszero(a[end]) && return a
	  pop!(a)
	  return trimlastzeros!(a)
	end



	function zplane(zs, ps; kwargs...)
		scatter(real.(zs), imag.(zs);
			  marker = (:black, :circle), label="Cero", kwargs...)
		scatter!( real.(ps), imag.(ps);
			marker = (:red, :xcross), label="Polo", kwargs...)
	  ts = range(0,stop=2pi;length=100)
	  plot!(cos.(ts), sin.(ts); aspect_ratio = 1, legend=false, kwargs...)
	end

	zplane(pr::DSP.PolynomialRatio; kwargs...) = 
		zplane(DSP.ZeroPoleGain(pr); kwargs...)
	
	DSP.filt(zpg::DSP.ZeroPoleGain, r...; kwargs...) = 
		filt(polynomialratio(zpg), r...; kwargs...)



	# Delta
	d(n) = n == 0 ? 1. : 0.

	# Escalón
	u(n) = n >= 0 ? 1. : 0.


end;

# ╔═╡ 8f1394ee-3e63-11eb-0093-e75468460dc5
md"# Trabajo Práctico Especial de Señales y Sistemas"

# ╔═╡ 7c04611e-3e61-11eb-0aa5-eb97132ace53
html"""
<h2 class="pm-node nj-subtitle">Identificación de canciones mediante huellas digitales acústicas</h2>
"""

# ╔═╡ adc46380-3e63-11eb-2422-5bfe1b5052ba
md"""
# Introducción

Imaginen la siguiente situación: se encuentran en un bar, tomando su trago favorito, inmersos en esa interesante conversación cuando de repente . . .  “pará, pará, ¿qué es eso que suena?”. Entre el ruido general alcanzan a distinguir esa canción que no escuchaban en tanto tiempo, que habían obsesivamente intentado encontrar pero que, a pesar de ello, nunca llegaron a dar siquiera con el nombre del intérprete. Su corazón se estremece . . .  Ahora que la tienen ahí, sonando nuevamente, la situación es bien diferente a la de hace algunos años: toman su teléfono celular, abren alguna de las aplicaciones de reconocimiento de audio que tienen y, en cuestión de segundos, aparece la información completa del tema, la banda, ¡e incluso se atreve a sugerir artistas similares!

Si el protagonista de esta ficción fuera estudiante de Señales y Sistemas, no debería extrañarnos que surjan en él varios interrogantes: ¿cómo hizo el programa para reconocer la canción escuchando sólo 10 segundos?, ¿cómo puede funcionar con el ruido de un bar de fondo? y ¿cómo pudee ser que además identifique tan rápido a seta banda que nadie conoce?

Resulta que todo este proceso se basa en reconocer *huellas digitales acústicas*, una técnica robusta y eficiente que puede entenderse en términos de Señales y Sistemas.

## Reconocimiento mediante huellas digitales acústicas

El reconocimiento de audio mediante huellas digitales acústicas es una técnica que busca identificar una pieza de audio, contrastando la información contenida en dicha pieza contra la almacenada en una base de datos de pistas conocidas. Esta técnica comienza a desarrollarse desde el año 2000, pero es durante la última década que tomé mayor impulso, siendo posible hoy en día encontrar una variedad de aplicaciones para smartphones capaces de identificar casi cualquier pieza de audio en segundos \[4\].

Existen varias formas dea bordar el problema del reconocimiento, pero todas persiguen la misma finalidad:

* Simplicidad computacional: el reconocimiento debe realizarse en forma rápida.
* Eficiencia en el uso de memoria y buena capacidad de discriminación: existen aproximadamente 100 millones de canciones grabadas \[1\].
* Robustez ante degradaciones: ruido de fondo, degradación por compresión digital del audio, ecualización por respuesta en frecuencia no plana del lugar y/o parlantes, etc.
* Granularidad: capacidad de reconocimiento utilizando sólo un segmento del audio.
* Alta tasa de aciertos, baja tasa de falsos positivos y de rechazos.

El procedimiento básico consiste en analizar el segmento de audio buscando extraer características particulares en el esapcio tiempo-frecuencia (es decir, en su espectrograma) que sirvan para identificarlo luego. Una característica podría ser, por ejemplo, la potencia que posean las diferentes bandas de frecuencias. Al conjunto de estas características se las denomina *huella digital acústica*.

![sys-tpe1.png](https://i.imgur.com/FoKz8Tw.png)

Este procedimiento de extracción de huellas se utiliza tanto para confeccionar la base de datos de canciones conocidas como para posteriormente reconocer segmentos de audio que se desean identificar consultando la base de datos. 

La elección de las características es la base del funcionamiento del método. Deben ser lo suficientemente específicas como para identificar a cada canción unívocamente pero a la vez ocupar poca memoria para permitir realizar la búsqueda en forma rápida y eficiente. También deberán ser relativamente inmunes ante ruidos y distorsiones del audio de manera de lograr un sistema de reconocimiento robusto. 

Diferentes características han sido propuestas con el fin de lograr estas especificaciones. Algunos ejemplos que se extraen del dominio tiempo-frecuencia son los MFCC \[4\] (Mel-Frequency Cepstrum Coefficients, comúnmente utilizados para la representación del habla), SFM \[1\] (Spectral Flatness Measure, una estimación de cuánto se aproxima una banda espectral a ser un tono o ruido), la ubicación de los picos de potencia del espectrograma (la base del algoritmo de Shazam \[9, 8\]), la energía de las bandas \[6\], etc. Además del espectrograma, otras técnicas se han utilizados como *wavelets* y algoritmos de visión \[2\], algoritmos de machine learning \[3\], y meta-características \[7\], útiles en grabaciones muy distorsionadas. 

En este trabajo desarrollaremos una implementación simple que utiliza como características el signo de la diferencia de energía entre bandas \[5\]. 

Cabe destacar que esta característica sólo sirve para reconocer las mismas versiones de las canciones que se almacenan en la base de datos: es decir, no está diseñada para reconocer versiones en vivo o interpretaciones diferentes a la original.

# Descripción del sistema de reconocimiento

El sistema de reconocimiento consta de dos bloques principales:

1. El algoritmo para extraer las huellas digitales acústicas.
2. La base de datos que almacena las huellas acústicas de las canciones conocidas y permite realizar búsquedas a partir de la huella acústica de una canción desconocida.

El algoritmo para extraer las huellas acústicas toma como entrada el archivo con el audio, lo pre-procesa, extrae las características, y entrega como resultado la huella acústica.


![sys-tpe2.png](https://cdn.nextjournal.com/data/QmRYiHJKoMuXyNybtpYDvExYoxBiwoA4dHQei9Quar5Sit?filename=sys-tpe2.png&content-type=image/png)

La base de datos guardará las huellas acústicas de las canciones conocidas. Tiene además un sistema de búsqueda (en inglés query) tal que al realizar un query con una huella acústica dada nos devuelve la canción – más probable – a la cual correseponde.

El esquema del sistema completo se presenta a continuación:

![sys-tpe3.png](https://cdn.nextjournal.com/data/QmYQPVLaSUHo9XBXL2m3HqMnbA2QAY9jWg7VqoFeby8N77?filename=sys-tpe3.png&content-type=image/png)

Observe que el algoritmo de extracción de huellas se utiliza tanto en la creación de la base de datos de canciones conocidas como para el reconocimiento de audios desconocidos.

## Algoritmo de extracción de huellas digitales acústicas

El algoritmo de extracción de huellas acústicas tiene como finalidad extraer características a partir del espectrograma del audio. Primeramente, se acondiciona el audio con un pre-procesamiento. Luego se realiza un espectrograma del audio pre-procesado y finalmente se extraen las características para cada intervalo de tiempo.

### 1. Pre-procesamiento del audio

Se comienza convirtiendo el audio, que suele ser estar en estéreo, a un audio monocanal, promediando los canales. 

Luego, se reduce su frecuencia de muestreo, aprovechando que son las frecuencias bajas las que contienen la información más útil para la extracción de características. En general, el audio está muestreado a 44100 Hz (calidad de CD), pero para este algoritmo de reconocimiento basta con tenerlo muestreado a 1/8 de su frecuencia original, 5512.5 Hz. Esto permite trabajar con menos muestras, aliviando la carga computacional. Para realizar los siguientes ejercicios, utilice el archivo `Pink.ogg`.

#### Ejercicio 1)

**Cargue la pista de audio. Verifique que la variable cargada es una matriz con 2 columnas, cada una correspondiendo a un canal de audio. Promedie los canales utilizando al función mean para tener 1 solo canal y grafique una porción de la señal resultante. Escuche el audio resultante para verificar el resultado.**
"""

# ╔═╡ a3bf22c4-3ea3-11eb-3d3d-adfdfc171c33
# La frecuencia de muestreo está fija a 44100 Hz
const sr = 44100

# ╔═╡ d132a762-3ea3-11eb-3494-692576a31f34
function loadaudio(fn)
	
    xraw = load(fn)
	
		# Chequeo que la frecuencia de muestreo sea `sr`
    @assert samplerate(xraw) == sr  samplerate(xraw) 

    return collect(xraw)
end

# ╔═╡ 28c5ed26-3e6b-11eb-1d44-01e209b92f00
# Wrapper para ver un vector de muestras de audio a 44100 Hz como un widget escuchable
sound(x) = SampleBuf(x, sr);

# ╔═╡ edf9177a-7015-11eb-0553-9512d2a99f22
pink = loadaudio("Pink.ogg")

# ╔═╡ 04aa1712-7016-11eb-10cc-99a4e40b5eef
pinkSound = sound(pink)

# ╔═╡ 2a80005a-7093-11eb-0d88-c5c2b7018ced
function meanAudio(x, dims = 2)
	xmean = mean(x, dims=dims)
	return [xmean[i] for i in 1:length(xmean)]
end;

# ╔═╡ 11802b5c-7016-11eb-1b8f-21b09c8c6137
pinkMono = meanAudio(pink)

# ╔═╡ 21359f50-7016-11eb-2198-ff85c9920a7b
pinkMeanSound = sound(pinkMono)

# ╔═╡ b9ad22ac-3e67-11eb-35e1-7f4579b64838
md"""
#### Ejercicio 2)

**Mediante la función `fft`, realice un análisis espectral de la señal mostrando la densidad espectral de potencia (PSD) en función de la frecuencia en Hz.  Verifique que la mayor PSD se concentra en las bajas frecuencias. Utilice escala logarítmica para mostrar esultado.**

**La PSD, ``S_{XX}``, la puede estimar como ``S_{X X}(\omega) = 1/T |X(\omega)|^2``, donde ``T`` es la duración del segmento temporal que usa para la estimación, o mejor aún, dividiendo al audio en segmentos, realizando la estimación anterior con cada uno, y finalmente promediándolas.**
"""

# ╔═╡ 50c411a8-7088-11eb-0a9d-2fe52da17a2a
function getpsd(x, rate)
	xrated = [x[i:i+rate-1] for i in 1:rate:length(x)]
	xratedfft = fft.(xrated)
	xratedfft = [(abs.(i)) .^2 ./ (sr * rate) for i in xratedfft]
	xratedfftmean = mean(xratedfft)
	xratedfftdb = log10.(xratedfftmean) .* 10
	return xratedfftdb[1:Int(round(length(xratedfftdb)/2))]
end;

# ╔═╡ 60c5245c-7088-11eb-1d2d-3f552a91ae07
psd = getpsd(pinkMono, 300)

# ╔═╡ 7d3d6748-7088-11eb-12c6-070b01843e4b
plot(
	range(0; stop=sr/2, length=length(psd)),
	psd;
	xlabel="Frecuencia (Hz)",
	ylabel = "PSD (dB)"
)

# ╔═╡ dc717174-7ad7-11eb-3f8a-97a269483324
md"""
Como se observa, la frecuencias mas bajas tienen mayor potencia.
"""

# ╔═╡ b60ae59e-3e67-11eb-123e-11c0cba7d09e
md"""
#### Ejercicio 3)

**Mediante la función `fft`, obtenga el espectro de una ventana rectangular y una de Hamming de igual duración y grafique su potencia en escala logarítmica en un mismo gráfico para compararlos. ¿Cuál es la resolución en frecuencia de cada ventana? Al realizar un espectrogama, ¿qué ventaja tiene utilizar la ventana de Hamming frente a la rectangular?**
"""

# ╔═╡ 8bce8660-794e-11eb-0748-c97c77e90cba
md"""
La ventana de hamming tiene lobulos laterales significativamente menores respecto a los lobulos laterales de la ventana rectangular, aunque el lobulo principal es aproximadamente el doble que la de la ventana rectangular.
Al tener lobulos laterales mas grandes, tendremos una mayor amplificacion en las frecuencias. Esto provoca que en el espectograma se muestren frecuencias que no correspondan.  
"""

# ╔═╡ 8d25b3bc-713d-11eb-19e5-e50944c500c0
begin
	samples = 20
	toDb(x) = 20log10(x) 
	
	vRec = ones(samples)
	vRecFFT = abs.(fft(padright(vRec,1000)))
	vRecFFT = fftshift(vRecFFT)
	vRecLog = toDb.(vRecFFT)
	
	vHamming = hamming(samples)
	vHammingFFT = abs.(fft(padright(vHamming,1000)))
	vHammingFFT = fftshift(vHammingFFT)
	vHammingLog = toDb.(vHammingFFT)
	
	plot(
		range(-pi; stop=pi, length=length(vRecLog)),
		vRecLog;
		xlabel="Ω",
		ylabel="dB",
		xlimit= (0, pi),
		ylimit= (-50, 30),
		label= "Rectangular",
		title="Espectro"
	)
	plot!(
		range(-pi; stop=pi, length=length(vHammingLog)),
		vHammingLog;
		xlabel="Ω",
		xlimit= (0, pi),
		label= "Hamming"
	)
end

# ╔═╡ 4aa62c04-7ae0-11eb-00e2-cd199c6cc45a
let
	plot(
		range(-sr/2; stop=sr/2, length=length(vRecLog)),
		vRecLog;
		xlabel="Frecuencia [Hz]",
		ylabel="dB",
		xlimit= (0, sr/2),
		ylimit= (-50, 30),
		label= "Rectangular",
		title= "Equivalencia en frecuencia"
	)
	plot!(
		range(-sr/2; stop=sr/2, length=length(vHammingLog)),
		vHammingLog;
		xlabel="Frecuencia [Hz]",
		xlimit= (0, sr/2),
		label= "Hamming"
	)
end

# ╔═╡ f16081a6-7abe-11eb-3ae5-c5af181dc83e
md"""
Observando los gráficos, podemos notar que el lóbulo principal de la ventana rectangular tiene una resolución de aproximadamente 2.5 KHz y la de la ventana de hamming 5KHz
"""

# ╔═╡ b2025250-3e67-11eb-39a2-73292bbf17c9
md"""
#### Ejercicio 4)
**Implemente un sistema para reducir la frecuencia de muestreo del audio, de 44100 Hz a 5512.5 Hz. Muestre un diagrama en bloques y justifique su diseño. Diseñe el filtro pasabajos mediante el método de ventaneo explicando y justificando las decisiones, de forma tal que su retardo sea menor a 1 ms. Graficar un diagrama de polos y ceros, respuesta en frecuencia (módulo y fase), respuesta al impulso, y retardo de grupo.**
"""

# ╔═╡ 8135a7d8-7aac-11eb-13fb-3977a7f4ceff
md"""
Para reducir la frecuencia de muestreo debemos pasar la señal por un filtro pasa bajos y luego decimar.

![filtro y decimacion.png](https://github.com/rinaldia92/sys-tpe/blob/master/filtro%20y%20decimacion.png?raw=true)

Como debemos reducirlo a un tamaño 8 veces menor utilizaremos una frecuencia de corte equivalente a sr/8 y N = 8.

Pero el filtro ideal no tiene soporte acotado, por lo cual debemos recortarlo utilizando el metodo de ventaneo.

![filtro soporte acotado.png](https://github.com/rinaldia92/sys-tpe/blob/master/filtro%20soporte%20acotado.png?raw=true)

Si bien ahora tenemos un soporte acotado, debemos realizar un desplazamiento a la derecha de la respuesta al impulso para que el sistema sea causal, pero esto nos genera un desfasaje lineal y un retardo de grupo.

Para calcular el retardo de grupo debo analizar la cantidad de muestras por segundo que hay en las canciones, dado que la frecuencia de muestreo es de 44100Hz sabemos que hay 44100 muestras por segundo (44.1 muestras cada milisegundo).
Al aplicar el metodo de ventaneo tenemos un sistema que no es causal, por lo cual debemos realizar un desplazamiento para convertirlo en causal y esto genera un delay equivalente a (N - 1) / 2, siendo N el orden del filtro. Dado que queremos un delay menor a 1 ms tenemos que cumplir con que (N - 1)/2 < 44.1 dando N < 89.2 (N natural). Se eligio N = 81
"""

# ╔═╡ a8892f14-70b4-11eb-26ff-11bb7d0f047e
begin
	H_ideal(Ω) = u(Ω + π / 8) - u(Ω - π / 8)
	h_ideal(n) = 1/8 * sinc(n / 8)
	h_rec = h_ideal.(-40:40) # pasandola por una ventana rect entre -40 y 40
	H_rec = fft(h_rec)
	hamming_vec = hamming(81)
	h_hamming = h_rec .* hamming_vec
	H_hamming = fft(h_hamming)
	subRate(x) = [x[8*i + 1] for i in 0:Int(round(length(x)/8))-1]
end;

# ╔═╡ 9087df4c-7ad6-11eb-1a95-db70cdcb1640
md"""
###### Mostramos el efecto de filtrar usando usando una ventana rectangular y luego submuestrear.
"""

# ╔═╡ 8c48e750-70a4-11eb-269e-939ec902c471
pinkMonoFiltRec = conv(pinkMono,h_rec)

# ╔═╡ 999f6532-70a4-11eb-0f3f-35b51be6777a
pinkSubratedRec = subRate(pinkMonoFiltRec)

# ╔═╡ 4c56bc6c-3ea4-11eb-01e7-7b26c1d054f0
SampleBuf(pinkSubratedRec, sr/8)

# ╔═╡ b79f83f0-7ad6-11eb-0da6-87aa0eb4dff8
md"""
###### Mostramos el efecto de filtrar usando usando una ventana de hamming y luego submuestrear.
"""

# ╔═╡ c89ab634-70ba-11eb-36f8-0526c0bd45d9
pinkMonoFiltHamming = conv(pinkMono,h_hamming)

# ╔═╡ 67a9c698-70c3-11eb-1e95-4566df3be24c
pinkSubratedHamming = subRate(pinkMonoFiltHamming)

# ╔═╡ 73ed9db0-70c3-11eb-0eee-71f0507bc207
SampleBuf(pinkSubratedHamming, sr/8)

# ╔═╡ becc85f6-7ad6-11eb-0dff-8d9c85c17965
md"""
###### Graficos correspondientes al filtro usando ventana de hamming.
"""

# ╔═╡ 841b1cb4-70c3-11eb-395b-0d24a544a58b
stem(h_hamming; xlabel="n", title= "Respuesta al impulso")

# ╔═╡ e6c7024a-7944-11eb-0b80-f14b78567e66
function groupDelay(x,y)
	z = [0.0 for i in 1:length(x)-1]
	for  i in 1:length(x)-1
		z[i] = (y[i+1] - y[i])/(x[i+1]-x[i])
	end
	return z .* -1
end

# ╔═╡ 58a04e18-7944-11eb-0912-6951ee323eef
begin
	plotRange = range(-pi; stop = pi, length = length(H_hamming))
	plotAbs = fftshift(abs.(H_hamming))
	plotPhase = unwrap(angle.(fftshift(H_hamming)))
	plotDelay = groupDelay(plotRange, plotPhase)
end
	

# ╔═╡ 4186458e-7ab6-11eb-157c-b3896fcddbb5
md"""
Aca observamos el modulo de la respuesta en frecuencia. Se nota que la forma es muy similar al filtro ideal.
"""

# ╔═╡ 8fc53be4-70c3-11eb-274c-5d427b8aece0
plot(
	plotRange,
	plotAbs;
	xlabel="Ω",
	xlims = (0,pi),
	legend=false,
	title="Modulo"
)

# ╔═╡ a2e71ffc-7ad3-11eb-12ee-e35a5d5b633a
md"""
En el gráfico de la fase, en sector entre 0 y pi/8, vemos que esta es lineal, lo cual se corresponde con el desfasaje que se le aplico al filtro ideal. El desplazamiento fue de 40 puntos a la derecha, por lo cual se espera que haya una pendiente de 40 lo cual se observa en el gráfico.
"""

# ╔═╡ 09d9b998-7603-11eb-2b46-a74dad9d0bc0
plot(
	plotRange,
	plotPhase;
	xlabel="Ω",
	xlims = (0,pi),
	legend=false,
	title= "Fase"
)

# ╔═╡ 27ff4ee8-7945-11eb-0887-493ca41e602c
plot(
	plotRange[1:length(plotRange)-1],
	plotDelay;
	xlabel="Ω",
	xlims = (0,pi),
	legend=false,
	title= "Retardo de grupo"
)

# ╔═╡ 026f30e6-7ab0-11eb-1e8e-c1d26da0bc1a
md"""
En el siguiente diagrama de Polos y Ceros observamos que tenemos polos en el origen y ceros sobre el origen y sobre el circulo unitario, a excepcion de la zona correspondiete a ω chicas. Esto tiene sentido ya que un filtro FIR si tiene polos los tiene en el origen y/o infinito y estamos realizando un filtro pasabajos por los cuales las bajas frecuencias deben pasar.

Por las características del filtro, se esperan tener 80 polos en cero y como la fase es lineal se esperan que los ceros sean 80, siendo pares conjugados. En este caso, tenemos 79 ya que el cero ubicado en el origen es doble.
"""

# ╔═╡ 047e53b2-70c4-11eb-090a-7bb200921492
let pr = polynomialratio(h_hamming, [1])
	zplane(getzeros(pr), getpoles(pr))
end

# ╔═╡ e7a268d6-7ab5-11eb-2cec-87ddc4950d37
md"""
Aca comparamos las PSD correspondientes a la señal original y la señal filtrada, como se observa, en las frecuencias bajas ambas señales tienen la misma potencia, pero llegado a un punto, la potencia de las frecuencias mas altas decae abruptamente en la señal filtrada con respecto a la original, lo cual tiene sentido ya que estabamos filtrando con un pasabajos.
"""

# ╔═╡ 014112be-70c6-11eb-1b88-dfa9b073f542
psdFiltered = getpsd(pinkMonoFiltHamming[1:length(pinkMono)], 300)

# ╔═╡ 5f21cdf2-796f-11eb-2a52-d33158a3f251
let
	plot(
		range(0; stop=sr/2, length=length(psd)),
		psd;
		xlabel="Frecuencia (Hz)",
		ylabel = "PSD (dB)"
	)
	plot!(
		range(0; stop=sr/2, length=length(psd)),
		psdFiltered;
		xlabel="Frecuencia (Hz)",
		ylabel = "PSD (dB)"
	)
end

# ╔═╡ af4f3da4-3e67-11eb-3cc6-3378e0c12667
md"""
#### Ejercicio 5)

**Realice un espectrograma de la señal original y la filtrada y verifique los efectos del filtrado. Indique el tipo y longitud de ventana utilizada.**
"""

# ╔═╡ 1a1f970a-7978-11eb-1104-b944d4d678e5
md"""
Para ambos casos se utilizar una ventana de hamming de 1024 puntos con un solapamiento del 90%.
"""

# ╔═╡ 6bd3e248-70c9-11eb-35d2-5930819fa625
specplot(pinkMono; fs=sr, overlap=0.9, window = hamming(1024))

# ╔═╡ 43b74df2-7977-11eb-1fa1-35e602e54789
md"""
Como se observa en el espectograma, vemos que se detectan frecuencias de hasta 15 KHz pero se nota que las frecuencias bajas predominan.
"""

# ╔═╡ 8e806244-70c9-11eb-26ce-3da3372d362e
specplot(pinkMonoFiltHamming; fs=sr, overlap=0.9, window = hamming(1024))

# ╔═╡ 356a700c-7978-11eb-12a8-41978621f6aa
md"""
En este caso, observamos que hay frecuencias de hasta 2700 Hz aproximadamente. Lo cual tiene sentido debido a que la filtramos por un pasabajos con un frecuencia de corte 8 veces mas baja que la frecuencia de muestreo.
"""

# ╔═╡ 982538c4-3e67-11eb-229e-dd2531a540d6
md"""
### Espectrograma

Los parámetros con los que se realice el espectrograma tienen una influencia directa sobre la extracción de las características. La elección de la ventana (tipo y longitud) es el parámetro principal, ya que debería lograr una resolución razonable tanto en tiempo como en frecuencia que permita distinguir las características espectrales y temporales con claridad.

#### Ejercicio 6)

**Realice un espectrograma de la señal sub-muestreada y muestre una imagen (con zoom) de alguna región donde se vean las características espectrales de la señal dada la ventana utilizada. Muestre y justifique cómo cambia la visualización de las características con 3 diferentes longitudes de ventana y comente qué longitud utilizará en el algoritmo. Realice las comparaciones con ventana rectangular y de Hamming.**
"""

# ╔═╡ fab10476-7abf-11eb-2d87-3f393cde1685
md"""
###### Espectogramas con ventana de hamming de 512, 1024 y 2048 puntos, con un solapamiento del 90%.
"""

# ╔═╡ 7fcc10bc-797b-11eb-321e-0736a8844282
specplot(pinkSubratedHamming; fs=sr/8, overlap=0.9, window = hamming(512), ylims=(0, 1e3), xlims=(3,5), title= "Hamming 512")

# ╔═╡ 39f5fc86-3ea4-11eb-37f3-25feb7d2aee6
specplot(pinkSubratedHamming; fs=sr/8, overlap=0.9, window = hamming(1024),ylims=(0, 1e3), xlims=(3,5), title= "Hamming 1024")

# ╔═╡ 798ef802-7149-11eb-2fd4-672a67b4fce1
specplot(pinkSubratedHamming; fs=sr/8, overlap=0.9, window = hamming(2048),ylims=(0, 1e3), xlims=(3,5), title= "Hamming 2048")

# ╔═╡ 2e530b82-7ac0-11eb-1e46-ed970d25b507
md"""
###### Espectogramas con ventana rectangular de 512, 1024 y 2048 puntos, con un solapamiento del 90%.
"""

# ╔═╡ 8b254136-797b-11eb-222c-572bbd32e034
specplot(pinkSubratedRec; fs=sr/8, overlap=0.9, window = ones(512), ylims=(0, 1e3), xlims=(3,5), title= "Rectangular 512")

# ╔═╡ 9337a2da-7942-11eb-1c8b-0d572ed4d139
specplot(pinkSubratedRec; fs=sr/8, overlap=0.9, window = ones(1024),ylims=(0, 1e3), xlims=(3,5), title= "Rectangular 1024")

# ╔═╡ 98469fe2-7942-11eb-1125-fb9dcda326c6
specplot(pinkSubratedRec; fs=sr/8, overlap=0.9, window = ones(2048), ylims=(0, 1e3), xlims=(3,5), title= "Rectangular 2048")

# ╔═╡ 3a9fd2e4-7ac0-11eb-3aa7-aba357ba7797
md"""
A simple vista, si comparamos los espectogramas correspondientes a la ventanas de 512 puntos podemos notar que en el espectograma que utiliza la ventana de hamming las frecuencias se encuentran bastante mas diferenciadas con respecto a la que utilizar la ventana rectangular.

Por otro lado, al utilizar una ventana mas grande se podrá visualizar de mejor manera el contenido en frecuencias perdiendo resolucion temporal, mientras que, cuanto mas chica es, mejora la resolucion temporal perdiendo resolucion espectral.

En el caso de la ventana de 2048 vemos que hay mucho solapamiento en las bandas de tiempo. Por lo cual nos conviene ver el espectograma con un rango mas grande en el tiempo.
"""

# ╔═╡ dfd795fa-7acd-11eb-398e-fbd55823326d
specplot(pinkSubratedHamming; fs=sr/8, overlap=0.9, window = hamming(2048),ylims=(0, 1e3), xlims=(1,9), title= "Hamming 2048")

# ╔═╡ 6a71b2d4-7ace-11eb-36a3-c9910278b1ff
md"""
Dado que a mayor ventana hay mayor resolución espectral y dado que en la comparación entre las ventanas de hamming y rectangular tenia mejor resolución la de hamming, se optó por utilizar una ventana de hamming de 2048 puntos para el algoritmo.
"""

# ╔═╡ 9309e284-3e67-11eb-1ab2-612f6c748c3b
md"""
### Extracción de características

La función provista `stft` permite obtener la matriz $s$ con las sucesivas DFT de los segmentos de la señal ventaneados -- lo que se grafica en un espectrograma.

Estas DFT nos devuelven muestras del espectro en frecuencias equiespaciadas en escala lineal – junto con los vectores de tiempos y frecuencias correspondientes a cada columna y fila. 

A partir de esto, necesitaremos la energía en bandas de frecuencia equiespaciadas en escala logarítmica (esto es similar a la escala psicoacústica Bark) debido a que así funciona el oído humano, que es un buen reconocedor de canciones – la relación entre frecuencias y posición de resonancia de la membrana basilar es aproximadamente logarítmica.

Debemos obtener una matriz de energías `E`, cuyas columnas, al igual que las de la matriz del espectrograma `s`, representan características espectrales del $n$-ésimo segmento de señal ventaneado (frame `n`). Para todo `n`, el elemento `E[m, n]` de nuestra nueva matriz debe contener la energía total de todos los coeficientes de `s[k, n]` asociados a frecuencias que caen dentro de la $k$-ésima banda de frecuencia.

#### Ejercicio 7)

**Divida al espectrograma en 21 bandas de frecuencias, equiespaciadas en escala logarítmica (es decir, el cociente entre frecuencias de bandas consecutivas debe ser constante). La frecuencia inferior de la primera banda debe ser 300 Hz y la superior de la última, 2 kHz.** 
	"""

# ╔═╡ e5cdf836-7add-11eb-0fc9-0ff239e2a7ef
md"""
###### Calculo de S
"""

# ╔═╡ 80498c84-74cd-11eb-3988-812f11a04268
function calcS(x; window = hamming(1000), ov = 0.5)
	nfft = length(window)
	overlap = Int(round(nfft*ov))
	S = stft(x; overlap=overlap, window=window, nfft=nfft)
	S = abs.(S).^2
	return S
end;

# ╔═╡ e0c8a45a-74cd-11eb-3502-e12a3b527240
begin
	window = hamming(2048)
	S = calcS(pinkSubratedHamming; window=window, ov=0.9)
	bands = exp.(range(log(300); stop=log(2000), length=22))
	freq = [i for i in range(0; stop = 2756.25, length = Int(round(2048/2)))]
end

# ╔═╡ 790fcb16-77cc-11eb-2d6e-fbd70848cd76
S

# ╔═╡ f3e07a98-7add-11eb-35bd-0d5d06f25857
md"""
###### Calculo de E
"""

# ╔═╡ c632d7ec-74a9-11eb-1133-3f64303d3fbd
function calcE(S, bands, freq)
	cantFilas, cantCol = size(S)
	E = [zeros(cantCol) for x in 1:21]
	j = 1
	while freq[j] < bands[1]
		j += 1
	end;
	for i in 1:length(bands)-1
		cant = 0
		lowFreq = bands[i]
		highFreq = bands[i+1]
		while freq[j] < lowFreq
			j += 1
		end;
		while j < length(freq) && freq[j] >= lowFreq && freq[j+1] < highFreq
			for k in 1:cantCol
				E[i][k] += S[cantFilas * (k - 1) + j]
			end;
			cant += 1
			j += 1
		end;
		E[i] ./ cant
	end;
	return E
end;

# ╔═╡ 07a7e9fe-74aa-11eb-08d0-e1c9faf07e9f
E = calcE(S, bands, freq)

# ╔═╡ 8deaf928-3e67-11eb-0327-31e0f74de814
md"""
Finalmente, las características se obtienen mediante una función que opera sobre E(m, n) según:


![sys-tpe4.png](https://i.imgur.com/dGoWXWE.png)

En palabras, $H[m, n]$ es 1 si la diferencia de energía entre las bandas $m$ y $m+1$ para el *frame* actual n es mayor a la del *frame* anterior $n-1$. Experimentalmente, se verficó que estas características son robustas ante varios tipos de procesamientos y distorsiones del audio\[5\].

#### Ejercicio 8)

**Implemente el algoritmo de extracción de características, calculando la huella. Debido al efecto borde, $H[m, n]$ debería resultar una matriz de 20 filas. Ejecute el algoritmo sobre un segmento de audio y muestre una imagen de la huella digital acústica obtenida.**
"""

# ╔═╡ 0385cce6-7ade-11eb-37c7-25a1ea2b936a
md"""
###### Calculo de H
"""

# ╔═╡ 6476e9fc-3ea4-11eb-3873-b765108f4bab
function calcH(E)
	f = length(E)-1
	c = length(E[1])
	H = Int.(zeros(f,c))
	for m in 1:length(E)-1
		for n in 2:length(E[m])
			dif1 = E[m+1][n] - E[m][n]
			dif2 = E[m+1][n-1] - E[m][n-1]
			H[f*(n-1)+m] = dif1 > dif2 ? 1 : 0
		end
	end
	return H
end;

# ╔═╡ e600a89c-719d-11eb-0887-2325103368e2
H = calcH(E)

# ╔═╡ 04acb5a6-719e-11eb-1aba-fdebe28d983f
plotHuella(h) = Gray.(h);

# ╔═╡ 2365113c-719e-11eb-3eb7-499d250a39a5
plotHuella(H)

# ╔═╡ 89743a62-3e67-11eb-209e-9b1f3cc84e34
md"""

## Confección de la base de datos

En principio, la base de datos simplemente debe

* Guardar, para cada *frame*
  * las características obtenidas
  * un identificador de la canción de la cual provino
  * el número de *frame* dentro de la canción.
* Permitir buscar el identificador de canción y el número de *frame* a partir de una característica, o informar si la característica buscada no se encuentra en la base de datos.

Así, cuando se desea identificar una música desconocida, se calculan las características de cada *frame*, se obtienen los identificadores de la base, y se opera con ellos de algún modo razonable para devolver la (o las) canciones más probables.

Ahora bien, más allá de este trabajo, es importante que este tipo de algoritmos escalen a bases de datos grandes y funcione rápido; para eso, es crucial que la la base de datos aproveche bien el espacio en memoria y permita una búsqueda eficiente. Por esto, se le proveen funciones que implementan una versión sencilla de una base de datos que cumple con estas características, que deberán ser entendidos.

Para confeccionar la base de datos, se utilizará la función *generar_DB*. Esta función requiere que el alumno haya definido otra función que, dado un archivo de audio como entrada, lo pre-procese, analice y extraiga su huella digital acústica en forma automática. 


#### Ejercicio 9)

**Encapsule su algoritmo de extracción de huellas acústicas en una función llamada `generar_DB` que reciba como entrada el path del archivo de audio y devuelva su huella digital acústica.**

```julia
\"""
	generar_huella(fname)

Devuelve la huella del audio en el archivo `fname`.

Tiene que estar muestreado a 44100 Hz.
\"""
function generar_huella(fname::String)
	#...

	return huella
end
```
"""

# ╔═╡ 2083e654-7708-11eb-028d-c1f6a5a8e8ef
function addNoiseToSong(x, noise)
	noiseArray = [rand() for i in 1:length(x)]
	alpha = var(x)/exp(noise/10)
	noiseArray = noiseArray ./ alpha
	noiseArrayVar = var(noiseArray)
	beta = sqrt(alpha/noiseArrayVar)
	noiseArray = noiseArray .* beta
	return x .+ noiseArray
end

# ╔═╡ 23dabfe4-7708-11eb-1960-536d21165c69
function trimSong(x, duration)
	pps = sr*duration
	len = length(x)
	cant = Int(len-pps)
	index = rand(1:cant)
	return x[index:Int(index+pps)]
end

# ╔═╡ 741304fe-3ea4-11eb-15e8-09908d98ecb3
function generar_huella(fname::String; trim = false, duration = 0, addNoise = false, noise = 0)
	song = loadaudio(fname)
	songMono = meanAudio(song)
	if trim
		songMono = trimSong(songMono, duration)
	end;
	if addNoise
		songMono = addNoiseToSong(songMono, noise)
	end;
	songMonoFiltHamming = conv(songMono,h_hamming)
	songSubratedHamming = subRate(songMonoFiltHamming)
	window = hamming(2048)
	S = calcS(songSubratedHamming; window = window, ov=0.9)
	bands = exp.(range(log(300); stop=log(2000), length=22))
	freq = [i for i in range(0; stop = 2756.25, length = Int(2048/2))]
	E = calcE(S, bands, freq)
	H = calcH(E)
	return H
end;

# ╔═╡ 855a7d2e-3e67-11eb-0f46-a5c786d5caf3
md"""

#### Ejercicio 10)

**Observe que la cantidad de elementos a guardar en la base de datos se incrementa conforme la longitud de las ventanas del espectrograma inicial disminuye, o el solapamiento entre ventanas se incrementa. Determine el solapamiento entre ventanas del espectrograma para obtener una densidad de aproximadamente 25 elementos por segundo y utilice este valor para el ejercicio siguiente.**
"""

# ╔═╡ 622dc82a-7983-11eb-2fbb-9fafc8055997
md"""
Para el calculo del overlap, nos basamos en la documentacion de la funcion specgram de matlab.
"""

# ╔═╡ bd10e958-7981-11eb-3b7e-0ba2cbbb344a
md"""
![calculo-overlap.png](https://github.com/rinaldia92/sys-tpe/blob/master/calculo-overlap.png?raw=true)
"""

# ╔═╡ 39062338-7549-11eb-0ccd-c1557057bed4
function calcOverlap(x, fps, window)
	ns = length(x)
	duration = 8*ns/sr
	nw = length(window)
	f = duration * fps
	no = (nw * f - ns)/(f-1)
	return no/nw # Para calcular el porcentaje
end

# ╔═╡ 98bbe152-797f-11eb-1d45-7d9702995a94
calcOverlap(pinkSubratedHamming, 25, hamming(2048))

# ╔═╡ 53691f34-7ad8-11eb-2fe4-db2ec2abe6d6
md"""
Se redondeará y se tomará como overlap 0.9
"""

# ╔═╡ 81717fc8-3e67-11eb-05fc-5bde46597f8a
md"""

#### Ejercicio 11)

**Ejecute la función `generar_DB` para confeccionar la base de datos completa de su lista de canciones. Utilice al menos 40 canciones para llenar la base de datos. Puede usar la lista de canciones provista, y/o usar una lista de canciones propia. (Recuerde verificar que la frecuencia de muestreo de sus canciones sea de 44100 Hz).**
"""

# ╔═╡ b91537ac-3ea4-11eb-14d6-d341c535d83e
begin
		# Cambiar por el nombre del subdirectorio donde estén sus canciones
	songsdir = "40songs"	
	songs = readdir(songsdir)
end

# ╔═╡ 73333e92-3e85-11eb-26b6-7f0309ef2ee9
@bind songpicked Select(songs)

# ╔═╡ feb5d512-3e85-11eb-0116-29e4d9539595
LocalResource(joinpath(songsdir, songpicked))

# ╔═╡ 0cf7ba9c-3e74-11eb-18e2-c38aa20f9e9a
# Código de `generar_DB`
begin
	
	"Guarda la información para identificar un frame"
	struct FrameID
		songID     :: UInt16
		frameindex :: UInt32 # índice dentro de la base de datos
	end
	
	"Genera la base de datos a partir de la lista de archivos de audio"
	function generar_DB(songs; dir="")
		
			# Inicializa la base, que es un vector de 2^20 vectores de FrameIDs
		db = [ copy(FrameID[]) for _ in 1:2^20]

			
		for songid in 1:length(songs) # Para cada canción
			
				# Genera la huella
			huella = generar_huella(joinpath(dir, songs[songid]))
			
				# La agrega a la base de datos
			for i in 1:size(huella, 2) # Para cada frame en la huella de la canción
					# Agregar el landmark (huella del frame) a la base
				addlandmark!(db, huella[:, i], songid, i)
			end
		end

		return db
	end
	
	
		
		# Agrega el landmark (huella de un frame) a la base de datos, asociado a la
		# canción con ID `songid` y número de frame `frameidx`.
	function addlandmark!(db, lm, songid, frameidx)
		fr = FrameID(songid, frameidx)
		idx = landmark_to_index(lm)

		push!(db[idx], fr)
	end
	
		# convertir un landmark -- vector de 20 bools -- en un entero 
		# con esa representación binaria
	landmark_to_index(h) = sum(2 .^ (0:19) .* h) + 1
end;

# ╔═╡ e9255b8c-3e74-11eb-2960-5d01b0c99b13
db = generar_DB(songs; dir=songsdir);

# ╔═╡ 7c7c1424-3e67-11eb-1da0-5dbad0171b20
md"""
# Test del algoritmo

En esta etapa final, pondremos a prueba la eficacia del algoritmo completo de reconocimiento. El procedimiento consistirá en obtener el porcentaje de aciertos del método al reconocer segmentos de audio, los cuales serán sometidos a distintas distorsiones.

Para esto se suministra la función *query_DB*, la cual recibe como entradas la base de datos y la huella acústica del segmento de audio a reconocer, y devuelve el ID del la canción que mejor coincide con la huella. Para entender cómo opera esta función, lea el Apéndice y el código.
"""

# ╔═╡ 415e32e6-3e76-11eb-17fa-23bd653fb975
function query_DB(db, huella::AbstractMatrix)
	
		# Cantidad de landmarks en la huella
    n = size(huella, 2)
		
		# Extraigo la información de todos los frames que encajaron con algún landmark
		# Vea que puede haber más frames que landmarks si hubo firmas repetidas
    frames = [ fr for lm in eachcol(huella) for fr in db[landmark_to_index(lm)]]

		# Función que toma un songID y devuelve un puntaje según cuántas veces aparece
		# como máximo en `frames` en algún intervalo de tiempo de `n`
    function score(sid)
			
			# Rescato los índices de frame asociados al songID pedido
        fids = [ fr.frameindex for fr in frames if fr.songID == sid ]
		
			# Armo un vector con deltas de peso 1 ubicadas en los índices de frames
        aux = zeros(Int, maximum(fids))
        aux[fids] .= 1
		
			# Convoluciono eso con una ventana de ancho `n` y devuelvo el máximo
        return maximum(conv(ones(Int, n), aux))
    end

		# Listo todos los songIDs candidatos, borrando duplicados
    sids = union(getfield.(frames, :songID))

		# Devuelvo el songID que maximice el `score`
    return sids[findmax(score.(sids))[2]]
end;

# ╔═╡ 76ce23dc-3e67-11eb-0be0-91b6781840fb
md"""

#### Ejercicio 12)

**Evalúe la tasa de aciertos del algoritmo identificando segmentos de duración $T$ con tiempo inicial elegido al azar de canciones elegidas al azar (vea la función `rand`). Las canciones deberán ser las mismas que utilizó para confeccionar la base de datos. Realice la evaluación para 50 segmentos con duración T  entre 5, 10 y 20 segundos cada vez (150 evaluaciones en total) obteniendo la tasa de aciertos en cada caso.**
"""

# ╔═╡ 2c813b30-7ade-11eb-1f1d-1b2dbcd3a4f3
md"""
Funcion para testear con distintas duraciones, con posibilidad de cambiar la cantidad y agregar ruido.
"""

# ╔═╡ 92d96440-7771-11eb-0324-4bc254a0deb0
function test(durations; cant = 50, noises = [nothing])
	matchs = zeros(length(durations), length(noises))
	addNoise = noises != [nothing]
	for nIndex in 1:length(noises)
		for dIndex in 1:length(durations)
			pos = 0
			for i in 1:cant
				songid = rand(1:length(songs))
				hSong = generar_huella(joinpath(songsdir, songs[songid]); trim = true, duration = durations[dIndex], addNoise = addNoise, noise = noises[nIndex])
				queryid = Int(query_DB(db, hSong))
				if songid == queryid
					pos += 1
				end;
			end;
			matchs[dIndex, nIndex] = pos / cant;
		end;
	end;
	return matchs
end;

# ╔═╡ 54f82180-7adb-11eb-0cac-1b127676fe0b
md"""
Resultados de testear 50 muestras de 5s, 10s y 20s sin ruido.
"""

# ╔═╡ ac315eac-7771-11eb-2cd9-914977facd0d
results = test([5,10,20])

# ╔═╡ cde76dda-7adb-11eb-0908-2fc75dee872b
md"""
Como se observa, se obtuvo una tasa de 100% de acierto.
"""

# ╔═╡ 7229577a-3e67-11eb-0c71-f383056175d1
md"""
#### Ejercicio 13)

**Repita el ejercicio 12 sumando ruido a los segmentos de audio. Utilice la función `randn` para generar las muestras de ruido. Evalúe tasa de aciertos para $SNR =0dB$, $10dB$ y $20dB$, mostrando sus resultados en una tabla para 9 combinaciones de longitud temporal y ruido. Nota: $SNR=10 log_{10}(P_X/P_N)$ donde $P_X$ es la potencia media de la señal sin ruido, y $P_N$ es la potencia media del ruido sumado a la señal. Para el cálculo de la potencia media puede utilizar la función `var`, que estima la varianza de una señal, ya que las señales de audio no deberían componente continua o valor medio.**
"""

# ╔═╡ 913c05a8-7adb-11eb-32ed-d37f27705386
md"""
Resultados de testear 50 muestras de 5s, 10s y 20s con 0dB, 10dB y 20dB de ruido.

Las filas corresponden a la duración y las columnas al ruido.
"""

# ╔═╡ cc43c0ac-7719-11eb-2086-8d17e6c2fb36
resultsWithNoise = test([5,10,20]; noises = [0,10,20])

# ╔═╡ c8c5da9e-7adb-11eb-0d68-4d222154385d
md"""
Como se observa, se obtuvo una alta tasa de acierto incluso al haberle agregado ruido.
"""

# ╔═╡ 6d76f2f2-3e67-11eb-04dc-0580a2072dda
md"""
#### Ejercicio 14) *(OPTATIVO)*

**Reproduzca y grabe un segmento de alguna de las canciones e intente reconocerla con su algoritmo. Puede utilizar dispositivos como el micrófono de su PC, su teléfono celular, una radio, etc. Comente los resultados obtenidos así como las condiciones de ruido de fondo y los dispositivos utilizados. (Recuerde que su función debe recibir audios a 44100 Hz.)**
"""

# ╔═╡ 685698fa-3e67-11eb-2698-937dd4801b5c
md"""
#### Ejercicio 15) *(OPTATIVO, SOLO PARA LOS MÁS ENTUSIASTAS)*

**Repita el ejercicio 12 para T=10 solamente, afectando los segmentos con otro tipos de distorsiones elegidas entre las siguientes:**

* **Saturación**
* **Cuantización**
* **Cambio de velocidad**
* **Ecualización**

**Si elige saturación y/o cuantización, muestre el efecto que ocasiona sobre el audio mediante un espectrograma. Justifique si estas distorsiones pueden considerarse o no sistemas LTI.**
"""

# ╔═╡ 62b03a84-3e67-11eb-3949-2dc573c7d956
md"""
#### Ejercicio 16) *(OPTATIVO; SÓLO PARA LOS MÁS OBSESIVOS)*

**Verifique cómo cambia la tasa de aciertos cuando:**

* **cambia el solapamiento**
* **incrementa el solapamiento sólo al momento de identificar pero no al armar la base de datos**
* **cambia la longitud de la ventana manteniendo la tasa de frames por segundo**
* **cambia el algoritmo de extracción de características por algún otro que se le ocurra**
"""

# ╔═╡ 562997ce-3e67-11eb-015a-d318429ed230
md"""
# Apéndice

## La estructura de la base de datos

La base de datos es un vector de $2^{20}$ elementos, cada uno de los cuales es un vector de tamaño dinámico. Para cada frame de la huella acústica, se genera un valor a guardar en estos vectores de tamaño dinámico. El índice en el que se guarda cada valor se obtiene pasando a decimal el número binario conformado por las características del frame en cuestión. Note que la cantidad de elementos del vector de la base, $2^{20}$ se debe a que cada frame de la huella acústica posee 20 elementos binarios.

Cada valor a guardar es un objeto de tipo `FrameID`, que consiste de un entero de 32 bits sin signo, que determina el número de frame, y un entero de 16 bits sin signo que guarda el índice de la canción. Observe que el máximo número de canciones distintas que se podrán almacenar en esta tabla es $2^{16}$. Obserbe también que múltiples FrameIDs de distintos frames pueden requerir ser guardados en la misma posición por tener la misma huella; por esto, cada elemento es un vector que puede crecer según se requiera.

Esta implementación puede ser optimizada pero es un balance razonable entre eficiencia y simpleza, suficiente para este trabajo.


## La búsqueda de la canción más probable

La búsqueda en la tabla *hash* se realiza partiendo de una huella acústica de entrada y devuelve el ID de la canción más probable a la cual corresponde como salida. 

Para cada *frame* de la huella se obtiene el índice del vector de la base de datos que debe consultarse, pasando la característica de cada *frame* de binario a decimal, del modo inverso que cuando se almacenan elementos en la base, y se extraen todos los elementos almacenados en esa posición. Cada uno de estos elementos posee el ID de la canción a la cual corresponde y el número del frame que le correspondía originalmente.

Una forma sencilla de decidir a qué canción corresponde la huella podría ser elegir, de todos los elementos que se extrajeron, el ID que aparezca el mayor número de veces.

 Un refinamiento al criterio anterior es, dado que para cada ID se dispone del número de *frame* del cual se extrajo, quedarse con la mayor cantidad de ID dentro de un intervalo de *frames* igual a la longitud de la huella que se está consultando. Esto es lo que está implementado en la función *query_DB*. 

Por ejemplo, supongamos que se realiza un query a la base de datos con una huella que posee una longitud de 5 frames y se obtienen los siguientes 7 elementos:


![sys-tpe7.png](https://cdn.nextjournal.com/data/QmWrqDXaxxMXmbQ56JmTwH8EFfAPYmkB6uGQRwtqgBCXtH?filename=sys-tpe7.png&content-type=image/png)

Bajo el primer criterio se declararía ganador al ID 7 dado que aparece mayor cantidad de veces, pero bajo el segundo criterio se decide por el ID 3, debido a que en el intervalo de 5 frames que tiene la huella de entrada el ID 7 aparece como máximo 2 veces y el ID 3 aparece 3 veces.

# Bibliografía

![sys-tpe8.png](https://cdn.nextjournal.com/data/QmcYW98tYgXiX9gm4dzpWvQEjBW8dbr5bnwHkdUwRWvZE2?filename=sys-tpe8.png&content-type=image/png)
"""

# ╔═╡ Cell order:
# ╠═09062294-3e5f-11eb-176f-dfcbf841f111
# ╟─8f1394ee-3e63-11eb-0093-e75468460dc5
# ╟─7c04611e-3e61-11eb-0aa5-eb97132ace53
# ╟─adc46380-3e63-11eb-2422-5bfe1b5052ba
# ╠═a3bf22c4-3ea3-11eb-3d3d-adfdfc171c33
# ╠═d132a762-3ea3-11eb-3494-692576a31f34
# ╠═28c5ed26-3e6b-11eb-1d44-01e209b92f00
# ╟─edf9177a-7015-11eb-0553-9512d2a99f22
# ╠═04aa1712-7016-11eb-10cc-99a4e40b5eef
# ╠═2a80005a-7093-11eb-0d88-c5c2b7018ced
# ╠═11802b5c-7016-11eb-1b8f-21b09c8c6137
# ╠═21359f50-7016-11eb-2198-ff85c9920a7b
# ╟─b9ad22ac-3e67-11eb-35e1-7f4579b64838
# ╠═50c411a8-7088-11eb-0a9d-2fe52da17a2a
# ╠═60c5245c-7088-11eb-1d2d-3f552a91ae07
# ╠═7d3d6748-7088-11eb-12c6-070b01843e4b
# ╟─dc717174-7ad7-11eb-3f8a-97a269483324
# ╟─b60ae59e-3e67-11eb-123e-11c0cba7d09e
# ╟─8bce8660-794e-11eb-0748-c97c77e90cba
# ╠═8d25b3bc-713d-11eb-19e5-e50944c500c0
# ╠═4aa62c04-7ae0-11eb-00e2-cd199c6cc45a
# ╟─f16081a6-7abe-11eb-3ae5-c5af181dc83e
# ╟─b2025250-3e67-11eb-39a2-73292bbf17c9
# ╠═8135a7d8-7aac-11eb-13fb-3977a7f4ceff
# ╠═a8892f14-70b4-11eb-26ff-11bb7d0f047e
# ╟─9087df4c-7ad6-11eb-1a95-db70cdcb1640
# ╟─8c48e750-70a4-11eb-269e-939ec902c471
# ╟─999f6532-70a4-11eb-0f3f-35b51be6777a
# ╟─4c56bc6c-3ea4-11eb-01e7-7b26c1d054f0
# ╟─b79f83f0-7ad6-11eb-0da6-87aa0eb4dff8
# ╟─c89ab634-70ba-11eb-36f8-0526c0bd45d9
# ╟─67a9c698-70c3-11eb-1e95-4566df3be24c
# ╟─73ed9db0-70c3-11eb-0eee-71f0507bc207
# ╟─becc85f6-7ad6-11eb-0dff-8d9c85c17965
# ╟─841b1cb4-70c3-11eb-395b-0d24a544a58b
# ╟─e6c7024a-7944-11eb-0b80-f14b78567e66
# ╟─58a04e18-7944-11eb-0912-6951ee323eef
# ╟─4186458e-7ab6-11eb-157c-b3896fcddbb5
# ╟─8fc53be4-70c3-11eb-274c-5d427b8aece0
# ╟─a2e71ffc-7ad3-11eb-12ee-e35a5d5b633a
# ╟─09d9b998-7603-11eb-2b46-a74dad9d0bc0
# ╟─27ff4ee8-7945-11eb-0887-493ca41e602c
# ╟─026f30e6-7ab0-11eb-1e8e-c1d26da0bc1a
# ╟─047e53b2-70c4-11eb-090a-7bb200921492
# ╟─e7a268d6-7ab5-11eb-2cec-87ddc4950d37
# ╟─014112be-70c6-11eb-1b88-dfa9b073f542
# ╟─5f21cdf2-796f-11eb-2a52-d33158a3f251
# ╟─af4f3da4-3e67-11eb-3cc6-3378e0c12667
# ╟─1a1f970a-7978-11eb-1104-b944d4d678e5
# ╟─6bd3e248-70c9-11eb-35d2-5930819fa625
# ╟─43b74df2-7977-11eb-1fa1-35e602e54789
# ╟─8e806244-70c9-11eb-26ce-3da3372d362e
# ╟─356a700c-7978-11eb-12a8-41978621f6aa
# ╟─982538c4-3e67-11eb-229e-dd2531a540d6
# ╟─fab10476-7abf-11eb-2d87-3f393cde1685
# ╟─7fcc10bc-797b-11eb-321e-0736a8844282
# ╟─39f5fc86-3ea4-11eb-37f3-25feb7d2aee6
# ╟─798ef802-7149-11eb-2fd4-672a67b4fce1
# ╟─2e530b82-7ac0-11eb-1e46-ed970d25b507
# ╟─8b254136-797b-11eb-222c-572bbd32e034
# ╟─9337a2da-7942-11eb-1c8b-0d572ed4d139
# ╟─98469fe2-7942-11eb-1125-fb9dcda326c6
# ╟─3a9fd2e4-7ac0-11eb-3aa7-aba357ba7797
# ╟─dfd795fa-7acd-11eb-398e-fbd55823326d
# ╟─6a71b2d4-7ace-11eb-36a3-c9910278b1ff
# ╟─9309e284-3e67-11eb-1ab2-612f6c748c3b
# ╟─e5cdf836-7add-11eb-0fc9-0ff239e2a7ef
# ╠═80498c84-74cd-11eb-3988-812f11a04268
# ╟─e0c8a45a-74cd-11eb-3502-e12a3b527240
# ╟─790fcb16-77cc-11eb-2d6e-fbd70848cd76
# ╟─f3e07a98-7add-11eb-35bd-0d5d06f25857
# ╠═c632d7ec-74a9-11eb-1133-3f64303d3fbd
# ╟─07a7e9fe-74aa-11eb-08d0-e1c9faf07e9f
# ╟─8deaf928-3e67-11eb-0327-31e0f74de814
# ╟─0385cce6-7ade-11eb-37c7-25a1ea2b936a
# ╠═6476e9fc-3ea4-11eb-3873-b765108f4bab
# ╟─e600a89c-719d-11eb-0887-2325103368e2
# ╟─04acb5a6-719e-11eb-1aba-fdebe28d983f
# ╟─2365113c-719e-11eb-3eb7-499d250a39a5
# ╟─89743a62-3e67-11eb-209e-9b1f3cc84e34
# ╟─2083e654-7708-11eb-028d-c1f6a5a8e8ef
# ╟─23dabfe4-7708-11eb-1960-536d21165c69
# ╠═741304fe-3ea4-11eb-15e8-09908d98ecb3
# ╟─855a7d2e-3e67-11eb-0f46-a5c786d5caf3
# ╟─622dc82a-7983-11eb-2fbb-9fafc8055997
# ╟─bd10e958-7981-11eb-3b7e-0ba2cbbb344a
# ╠═39062338-7549-11eb-0ccd-c1557057bed4
# ╟─98bbe152-797f-11eb-1d45-7d9702995a94
# ╟─53691f34-7ad8-11eb-2fe4-db2ec2abe6d6
# ╟─81717fc8-3e67-11eb-05fc-5bde46597f8a
# ╟─b91537ac-3ea4-11eb-14d6-d341c535d83e
# ╟─73333e92-3e85-11eb-26b6-7f0309ef2ee9
# ╟─feb5d512-3e85-11eb-0116-29e4d9539595
# ╟─0cf7ba9c-3e74-11eb-18e2-c38aa20f9e9a
# ╠═e9255b8c-3e74-11eb-2960-5d01b0c99b13
# ╟─7c7c1424-3e67-11eb-1da0-5dbad0171b20
# ╟─415e32e6-3e76-11eb-17fa-23bd653fb975
# ╟─76ce23dc-3e67-11eb-0be0-91b6781840fb
# ╟─2c813b30-7ade-11eb-1f1d-1b2dbcd3a4f3
# ╠═92d96440-7771-11eb-0324-4bc254a0deb0
# ╟─54f82180-7adb-11eb-0cac-1b127676fe0b
# ╠═ac315eac-7771-11eb-2cd9-914977facd0d
# ╟─cde76dda-7adb-11eb-0908-2fc75dee872b
# ╟─7229577a-3e67-11eb-0c71-f383056175d1
# ╟─913c05a8-7adb-11eb-32ed-d37f27705386
# ╠═cc43c0ac-7719-11eb-2086-8d17e6c2fb36
# ╟─c8c5da9e-7adb-11eb-0d68-4d222154385d
# ╟─6d76f2f2-3e67-11eb-04dc-0580a2072dda
# ╟─685698fa-3e67-11eb-2698-937dd4801b5c
# ╟─62b03a84-3e67-11eb-3949-2dc573c7d956
# ╟─562997ce-3e67-11eb-015a-d318429ed230
