<h1 align="center">fortran-ballistics-sim</h1>

<p align="center"><b>Simulasi Numerik Gerak Peluru Dua Dimensi dengan Hambatan Udara Kuadratik<br>Menggunakan Metode Runge-Kutta Orde Empat dalam Fortran 2008</b></p>

<p align="center">NuRichter</p>

**_Abstrak_—** Repositori ini menyajikan simulasi gerak peluru dua dimensi di bawah
gravitasi seragam dan gaya hambat aerodinamis kuadratik. Sistem persamaan diferensial
biasa nonlinear orde satu diintegrasikan dengan metode Runge-Kutta orde empat (RK4)
berlangkah tetap dalam aritmetika presisi ganda IEEE 754, dilengkapi lokasi kejadian
berbasis bisection untuk menentukan titik puncak dan titik tumbukan secara tepat.
Implementasi diverifikasi terhadap solusi analitik vakum pada lima sudut elevasi dengan
galat relatif maksimum 2,4 × 10⁻¹³, dan studi konvergensi menghasilkan orde empiris
4,02 yang sesuai dengan orde teoretis metode. Untuk bola berdiameter 74 mm dan bermassa
145 g yang diluncurkan pada 50 m/s dengan sudut 45°, hambatan udara mengurangi jarak
tempuh dari 254,93 m menjadi 105,96 m (−58,4 %) dan menggeser sudut elevasi optimum
dari 45° ke 38,5°.

**_Kata Kunci_—** gerak peluru, hambatan udara kuadratik, Runge-Kutta orde empat,
lokasi kejadian, Fortran 2008, verifikasi numerik.

---

## I. Pendahuluan

Gerak peluru tanpa hambatan udara memiliki solusi tertutup berupa parabola [1].
Penambahan gaya hambat yang sebanding dengan kuadrat laju membuat komponen horizontal
dan vertikal saling bergantung secara nonlinear, sehingga solusi analitik umum tidak
tersedia dan lintasan harus dihitung secara numerik [1].

Kontribusi repositori ini adalah: (i) implementasi RK4 presisi ganda dalam Fortran 2008
standar tanpa pustaka eksternal; (ii) lokasi kejadian yang membuat besaran ringkasan
(jarak, tinggi maksimum, waktu tempuh) tidak dibatasi resolusi langkah waktu;
(iii) protokol verifikasi otomatis terhadap solusi analitik beserta studi konvergensi;
dan (iv) skrip visualisasi Python untuk keluaran lintasan.

## II. Metode

### A. Model Fisis

Peluru dimodelkan sebagai massa titik $m$ yang bergerak pada bidang vertikal $x$–$y$.
Gaya hambat aerodinamis

$$
F_d = \tfrac{1}{2}\,\rho\,C_d\,A\,v^2, \qquad A = \tfrac{\pi d^2}{4}, \qquad v = \sqrt{v_x^2 + v_y^2}
\tag{1}
$$

bekerja berlawanan arah vektor kecepatan, dengan $\rho$ massa jenis udara, $C_d$ koefisien
hambat, dan $d$ diameter peluru. Didefinisikan konstanta hambat per satuan massa

$$
k = \frac{\rho\,C_d\,A}{2m} \quad [\mathrm{m^{-1}}].
\tag{2}
$$

Hukum Newton kedua menghasilkan sistem orde satu untuk vektor keadaan
$\mathbf{s} = [x,\ y,\ v_x,\ v_y]^\mathsf{T}$:

$$
\dot{x} = v_x, \qquad
\dot{y} = v_y, \qquad
\dot{v}_x = -k\,v\,v_x, \qquad
\dot{v}_y = -g - k\,v\,v_y,
\tag{3}
$$

dengan syarat awal $\mathbf{s}(0) = [0,\ y_0,\ v_0\cos\theta,\ v_0\sin\theta]^\mathsf{T}$.
Integrasi dihentikan pada kejadian $y(t) = 0$ untuk $t > 0$. Model mengasumsikan
$\rho$, $C_d$, dan $g$ konstan, tanpa angin, tanpa gaya Magnus, dan permukaan tanah datar.
Nilai $\rho$ mengikuti atmosfer standar di permukaan laut [2], dan $C_d = 0{,}47$ adalah
nilai tipikal bola mulus pada rezim subkritis [3]. Persamaan (3) menyiratkan laju
terminal $v_t = \sqrt{g/k}$.

### B. Integrasi Numerik

Sistem (3), ditulis $\dot{\mathbf{s}} = \mathbf{f}(\mathbf{s})$, diintegrasikan dengan
metode Runge-Kutta orde empat klasik [4], [5] berlangkah tetap $h$:

$$
\begin{aligned}
\mathbf{K}_1 &= \mathbf{f}(\mathbf{s}_n), &
\mathbf{K}_2 &= \mathbf{f}\big(\mathbf{s}_n + \tfrac{h}{2}\mathbf{K}_1\big), \\
\mathbf{K}_3 &= \mathbf{f}\big(\mathbf{s}_n + \tfrac{h}{2}\mathbf{K}_2\big), &
\mathbf{K}_4 &= \mathbf{f}\big(\mathbf{s}_n + h\,\mathbf{K}_3\big),
\end{aligned}
\tag{4}
$$

$$
\mathbf{s}_{n+1} = \Phi_h(\mathbf{s}_n) = \mathbf{s}_n + \tfrac{h}{6}\big(\mathbf{K}_1 + 2\mathbf{K}_2 + 2\mathbf{K}_3 + \mathbf{K}_4\big).
\tag{5}
$$

Galat lokal metode ini $\mathcal{O}(h^5)$ dan galat global $\mathcal{O}(h^4)$ [4].
RK4 dipilih karena eksplisit, cukup akurat untuk sistem tak-kaku ini dengan empat
evaluasi $\mathbf{f}$ per langkah, dan tidak memerlukan iterasi implisit seperti skema
Verlet ketika gaya bergantung pada kecepatan. Untuk $C_d = 0$, solusi eksak (3) berupa
polinomial derajat dua dalam $t$ yang diintegrasikan RK4 tanpa galat pemotongan,
sehingga perbandingan dengan solusi vakum menguji kebenaran implementasi, bukan
akurasi diskretisasi. Waktu dihitung sebagai $t_n = n\,h$ untuk mencegah akumulasi galat
pembulatan.

### C. Lokasi Kejadian

Jika $y_n > 0$ dan $y_{n+1} \le 0$, dicari langkah parsial $\eta^\ast$ yang memenuhi

$$
\big[\Phi_\eta(\mathbf{s}_n)\big]_y = 0, \qquad \eta \in (0,\ h],
\tag{6}
$$

dengan bisection hingga interval tidak dapat dibagi lagi dalam presisi mesin [6].
Prosedur yang sama terhadap komponen $v_y$ memberikan titik puncak. Karena $\Phi_\eta$
adalah langkah RK4 itu sendiri, keadaan pada kejadian memiliki akurasi yang sama dengan
integrator, sehingga $R$, $H$, dan $T$ tidak dibatasi resolusi $h$.

### D. Implementasi

Kode ditulis dalam Fortran 2008 [7] dengan `implicit none` dan tipe `real64` dari modul
intrinsik `iso_fortran_env`, yaitu format biner64 IEEE 754 dengan
$\varepsilon \approx 2{,}2 \times 10^{-16}$ [8]. Modul `ballistics_core` memuat ruas kanan
(3), langkah RK4 (4)–(5), lokasi kejadian (6), dan solusi analitik vakum; program
`ballistics` menangani simulasi utama, verifikasi, dan studi numerik. Kompilasi dengan
gfortran 16.2.0 memakai opsi `-Wall -Wextra -Wimplicit-interface -fcheck=all -O2`
menghasilkan nol peringatan. Parameter bawaan tercantum pada Tabel I.

<p align="center">TABEL I<br>PARAMETER SIMULASI BAWAAN</p>

| Simbol | Parameter | Nilai | Satuan | Nama di kode |
|:---:|---|---:|:---:|---|
| $v_0$ | kecepatan awal | 50 | m/s | `v0` |
| $\theta$ | sudut elevasi | 45 | ° | `theta_deg` |
| $y_0$ | ketinggian awal | 0 | m | `y0` |
| $m$ | massa | 0,145 | kg | `mass` |
| $d$ | diameter | 0,074 | m | `diameter` |
| $C_d$ | koefisien hambat | 0,47 | – | `cd` |
| $\rho$ | massa jenis udara | 1,225 | kg/m³ | `rho` |
| $g$ | percepatan gravitasi | 9,80665 | m/s² | `g` |
| $h$ | langkah waktu | 1 × 10⁻³ | s | `dt` |
| $k$ | konstanta hambat (2) | 8,539 × 10⁻³ | m⁻¹ | turunan |

### E. Protokol Verifikasi

Tiga uji dijalankan otomatis oleh `./ballistics --verify`, dan program keluar dengan
kode 1 bila salah satunya gagal:

- **V1 (solusi vakum).** Dengan $C_d = 0$ dan $\theta \in \lbrace 15^\circ, 30^\circ, 45^\circ, 60^\circ, 75^\circ \rbrace$,
  hasil simulasi dibandingkan dengan

  $$
  R = \frac{v_0^2 \sin 2\theta}{g}, \qquad
  H_{\max} = \frac{v_0^2 \sin^2\theta}{2g}, \qquad
  T = \frac{2 v_0 \sin\theta}{g}.
  \tag{7}
  $$

  Kriteria lulus: galat relatif ketiga besaran $\le 10^{-10}$.
- **V2 (arah efek hambat).** Untuk $C_d = 0{,}47$, harus berlaku $R < R_\text{vakum}$ dan $H < H_\text{vakum}$.
- **V3 (konvergensi dua titik).** Selisih relatif $R$ antara $h = 10^{-2}$ s dan $h = 10^{-3}$ s harus $\le 10^{-6}$.

Studi konvergensi lengkap (`./ballistics --study`) menghitung galat relatif $e(h)$
terhadap solusi acuan $h_\text{ref} = 10^{-4}$ s untuk $h = 0{,}2/2^j$, $j = 0, \dots, 4$,
dengan orde empiris

$$
p = \log_2 \frac{e(h)}{e(h/2)}.
\tag{8}
$$

Mode yang sama juga mencari sudut elevasi yang memaksimalkan $R$ melalui sapuan
$\theta \in [25^\circ, 45^\circ]$ dengan resolusi $0{,}05^\circ$.

## III. Hasil

### A. Verifikasi terhadap Solusi Analitik

Tabel II menunjukkan hasil V1. Galat relatif berada pada orde 10⁻¹⁴–10⁻¹³, setara
dengan akumulasi pembulatan biner64 selama ribuan langkah. Hasil ini konsisten dengan
sifat RK4 yang eksak untuk percepatan konstan. Uji V2 dan V3 juga lulus, dengan
selisih relatif V3 sebesar 7,2 × 10⁻¹².

<p align="center">TABEL II<br>PERBANDINGAN SIMULASI (C<sub>d</sub> = 0) DENGAN SOLUSI ANALITIK VAKUM, v<sub>0</sub> = 50 m/s, h = 1 ms</p>

| $\theta$ | $R_\text{sim}$ (m) | $R_\text{analitik}$ (m) | $\lvert\Delta R\rvert / R$ | $\lvert\Delta H\rvert / H$ | $\lvert\Delta T\rvert / T$ |
|---:|---:|---:|---:|---:|---:|
| 15° | 127,464526622 | 127,464526622 | 4,6 × 10⁻¹⁴ | 1,9 × 10⁻¹⁴ | 1,6 × 10⁻¹⁴ |
| 30° | 220,775036272 | 220,775036272 | 3,9 × 10⁻¹⁴ | 8,8 × 10⁻¹⁴ | 4,7 × 10⁻¹⁴ |
| 45° | 254,929053244 | 254,929053244 | 1,9 × 10⁻¹³ | 1,2 × 10⁻¹³ | 8,1 × 10⁻¹⁴ |
| 60° | 220,775036272 | 220,775036272 | 5,5 × 10⁻¹⁴ | 1,4 × 10⁻¹³ | 9,5 × 10⁻¹⁴ |
| 75° | 127,464526622 | 127,464526622 | 2,4 × 10⁻¹³ | 1,4 × 10⁻¹³ | 1,0 × 10⁻¹³ |

### B. Konvergensi

Tabel III memuat galat terhadap solusi acuan untuk kasus berhambat. Orde empiris (8)
mendekati 4 dari atas, turun dari 4,13 ke 4,02 seiring mengecilnya $h$, sesuai orde
teoretis RK4. Ekstrapolasi $e \propto h^4$ memberikan galat diskretisasi sekitar
10⁻¹⁵ pada $h = 1$ ms, yaitu di bawah lantai pembulatan. Dengan demikian, ketidakpastian
pada hasil Bagian III-C didominasi oleh asumsi model, bukan oleh metode numerik.

<p align="center">TABEL III<br>GALAT RELATIF TERHADAP SOLUSI ACUAN h<sub>ref</sub> = 10⁻⁴ s DAN ORDE EMPIRIS, C<sub>d</sub> = 0,47</p>

| $h$ (s) | Langkah | $e_R$ | $e_H$ | $e_T$ | $p_R$ | $p_H$ | $p_T$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,2000 | 28 | 1,36 × 10⁻⁶ | 3,10 × 10⁻⁶ | 1,12 × 10⁻⁶ | – | – | – |
| 0,1000 | 55 | 7,77 × 10⁻⁸ | 1,82 × 10⁻⁷ | 6,45 × 10⁻⁸ | 4,13 | 4,10 | 4,12 |
| 0,0500 | 109 | 4,65 × 10⁻⁹ | 1,10 × 10⁻⁸ | 3,86 × 10⁻⁹ | 4,06 | 4,05 | 4,06 |
| 0,0250 | 218 | 2,84 × 10⁻¹⁰ | 6,74 × 10⁻¹⁰ | 2,36 × 10⁻¹⁰ | 4,03 | 4,03 | 4,03 |
| 0,0125 | 436 | 1,75 × 10⁻¹¹ | 4,18 × 10⁻¹¹ | 1,46 × 10⁻¹¹ | 4,02 | 4,01 | 4,02 |

### C. Lintasan dengan Hambatan Udara

Gbr. 1 membandingkan lintasan berhambat dengan lintasan vakum dari syarat awal yang
sama, dan Tabel IV merangkum besaran utamanya. Hambatan udara mengurangi jarak 58,4 %
dan tinggi maksimum 42,2 %. Lintasan menjadi tidak simetris: puncak terletak pada
58,0 % jarak tempuh, dan fase turun (2,976 s) lebih lama daripada fase naik (2,464 s).
Pada saat peluncuran, gaya hambat bernilai 2,18 kali berat peluru karena $v_0$ melampaui
laju terminal $v_t = 33{,}89$ m/s. Laju tumbukan 25,06 m/s masih di bawah $v_t$, sehingga
peluru belum mencapai kecepatan terminal saat menyentuh tanah. Sapuan sudut memberikan
jarak maksimum 107,918 m pada $\theta^\ast = 38{,}50^\circ$, yaitu 1,85 % lebih jauh
daripada pada 45°.

<p align="center">
  <img src="results/trajectory.png" width="760" alt="Lintasan peluru dengan hambatan udara dan lintasan vakum pembanding">
</p>
<p align="center"><sub>Gbr. 1. Lintasan peluru dengan hambatan udara (garis biru, RK4, h = 1 ms) dan lintasan vakum dari syarat awal yang sama (garis putus-putus). Titik menandai puncak dan titik tumbukan.</sub></p>

<p align="center">TABEL IV<br>BESARAN LINTASAN DENGAN DAN TANPA HAMBATAN UDARA (PARAMETER TABEL I)</p>

| Besaran | Simbol | Dengan hambatan | Vakum | Perubahan |
|---|:---:|---:|---:|---:|
| Jarak tempuh | $R$ | 105,958 m | 254,929 m | −58,4 % |
| Tinggi maksimum | $H$ | 36,861 m | 63,732 m | −42,2 % |
| Posisi horizontal puncak | $x_H$ | 61,434 m | 127,465 m | −51,8 % |
| Waktu ke puncak | $t_H$ | 2,464 s | 3,605 s | −31,7 % |
| Waktu tempuh | $T$ | 5,440 s | 7,210 s | −24,6 % |
| Laju tumbukan | $v_\text{imp}$ | 25,055 m/s | 50,000 m/s | −49,9 % |

## IV. Pembahasan

Dengan viskositas udara $\mu = 1{,}81 \times 10^{-5}$ Pa·s, bilangan Reynolds
$\mathrm{Re} = \rho v d/\mu$ selama penerbangan berada pada rentang
$9{,}1 \times 10^4$–$2{,}50 \times 10^5$. Batas bawah terjadi pada laju minimum
18,10 m/s ($t = 3{,}01$ s, sesudah puncak) dan batas atas pada saat peluncuran.
Rentang ini mencakup rezim subkritis hingga batas krisis hambat bola mulus
($\mathrm{Re} \approx 2$–$3 \times 10^5$), tempat $C_d$ dapat turun tajam [3]. Oleh
karena itu, asumsi $C_d$ konstan merupakan sumber ketidakpastian model terbesar, jauh
melampaui galat numerik pada Bagian III-B. Pengembangan yang wajar mencakup $C_d(\mathrm{Re})$,
$\rho(y)$ menurut atmosfer standar [2], angin, dan gaya Magnus untuk peluru yang berputar.

## V. Kesimpulan

Integrator RK4 presisi ganda dengan lokasi kejadian telah diimplementasikan dan
diverifikasi. Implementasi mereproduksi solusi vakum hingga batas pembulatan biner64,
menunjukkan orde konvergensi empiris 4,02, dan untuk parameter Tabel I memperlihatkan
pengurangan jarak 58,4 % serta pergeseran sudut optimum ke 38,5° akibat hambatan udara
kuadratik.

## Referensi

[1] J. R. Taylor, *Classical Mechanics*. Sausalito, CA, USA: University Science Books, 2005, ch. 2.

[2] *Standard Atmosphere*, ISO 2533:1975, International Organization for Standardization, Geneva, Switzerland, 1975.

[3] F. M. White, *Fluid Mechanics*, 7th ed. New York, NY, USA: McGraw-Hill, 2011, ch. 7.

[4] J. C. Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd ed. Chichester, U.K.: Wiley, 2016.

[5] E. Hairer, S. P. Nørsett, and G. Wanner, *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed. Berlin, Germany: Springer-Verlag, 1993.

[6] L. F. Shampine and S. Thompson, "Event location for ordinary differential equations," *Comput. Math. Appl.*, vol. 39, no. 5–6, pp. 43–54, 2000.

[7] *Information Technology—Programming Languages—Fortran—Part 1: Base Language*, ISO/IEC 1539-1:2010, 2010.

[8] *IEEE Standard for Floating-Point Arithmetic*, IEEE Std 754-2019, 2019.

---

## Lampiran A. Reproduksi

### Struktur Direktori

```
fortran-ballistics-sim/
├── src/
│   └── ballistics.f90        # modul ballistics_core + program utama
├── scripts/
│   └── plot_trajectory.py    # visualisasi results/trajectory.csv (Gbr. 1)
├── results/
│   ├── .gitkeep
│   └── trajectory.png        # Gbr. 1 (CSV tidak di-commit)
├── .gitignore
└── README.md
```

### Prasyarat

| Perangkat | Versi | Instalasi |
|---|---|---|
| gfortran | ≥ 9 | Linux: `sudo apt install gfortran` · macOS: `brew install gcc` · Windows: `winget install BrechtSanders.WinLibs.POSIX.UCRT` |
| Python | ≥ 3.9 | `pip install matplotlib pandas numpy` |

### Kompilasi dan Eksekusi

Jalankan dari akar proyek, karena program menulis ke `results/trajectory.csv`.
Di Windows, ganti `./ballistics` dengan `.\ballistics.exe`.

```bash
gfortran -Wall -Wextra -Wimplicit-interface -fcheck=all -O2 src/ballistics.f90 -o ballistics

./ballistics              # simulasi utama, ringkasan + results/trajectory.csv  (Tabel IV)
./ballistics --verify     # uji V1–V3                                          (Tabel II)
./ballistics --study      # konvergensi dan sudut optimum                       (Tabel III, III-C)
python scripts/plot_trajectory.py   # results/trajectory.png                   (Gbr. 1)
```

Opsi skrip plot: `--input`, `--output`, `--dpi`, dan `--no-vacuum`.

### Format Keluaran

`results/trajectory.csv` berisi satu baris per langkah waktu dengan 17 digit
signifikan, dan baris terakhir adalah titik tumbukan ($y = 0$).

| Kolom | Arti | Satuan |
|---|---|---|
| `t` | waktu | s |
| `x`, `y` | posisi horizontal dan vertikal | m |
| `vx`, `vy` | komponen kecepatan | m/s |

### Mengubah Parameter

Semua parameter pada Tabel I berada pada tipe `sim_params` di bagian atas
`src/ballistics.f90` (blok `PARAMETER AWAL`). Ubah nilainya lalu kompilasi ulang.
Jika `g` diubah, samakan pula konstanta `G` di `scripts/plot_trajectory.py` agar
lintasan vakum pada Gbr. 1 tetap benar. Rumus pembanding (7) hanya berlaku untuk $y_0 = 0$.
