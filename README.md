# fortran-ballistics-sim

Simulasi gerak peluru dua dimensi dengan hambatan udara kuadratik, ditulis dalam
Fortran modern (Fortran 2008) dengan presisi ganda. Persamaan gerak diintegrasikan
memakai metode Runge-Kutta orde 4 (RK4), hasil lintasan diekspor ke CSV, lalu
divisualisasikan dengan Python.

![Lintasan peluru dengan dan tanpa hambatan udara](results/trajectory.png)

## Model fisika

Peluru dianggap massa titik bermassa $m$ yang bergerak pada bidang $x$–$y$ di bawah
gravitasi konstan $g$ dan gaya hambat udara kuadratik:

$$
F_d = \tfrac{1}{2}\,\rho\,C_d\,A\,v^2, \qquad A = \tfrac{\pi d^2}{4}, \qquad v = \sqrt{v_x^2 + v_y^2}
$$

Gaya hambat selalu berlawanan arah dengan kecepatan, sehingga
$\vec F_d = -\tfrac{1}{2}\rho C_d A\,v\,\vec v$. Dengan $k = \dfrac{\rho\,C_d\,A}{2m}$
(satuan m⁻¹), sistem persamaan diferensial biasa orde satu menjadi:

$$
\frac{dx}{dt} = v_x, \qquad
\frac{dy}{dt} = v_y, \qquad
\frac{dv_x}{dt} = -k\,v\,v_x, \qquad
\frac{dv_y}{dt} = -g - k\,v\,v_y
$$

Kondisi awal: $x_0 = 0$, $y_0$ (bawaan 0), $v_{x0} = v_0\cos\theta$, $v_{y0} = v_0\sin\theta$.
Karena suku $v\,v_x$ dan $v\,v_y$ membuat gerak horizontal dan vertikal saling
bergantung secara nonlinear, sistem ini tidak memiliki solusi analitik tertutup dan
harus diselesaikan secara numerik.

Asumsi: tidak ada angin, $\rho$ dan $C_d$ konstan (tidak bergantung ketinggian atau
bilangan Reynolds/Mach), tanpa rotasi peluru (efek Magnus), dan bumi datar.

## Metode numerik

- **RK4 langkah tetap.** Galat global $\mathcal{O}(\Delta t^4)$ dengan empat evaluasi
  fungsi per langkah. Untuk $C_d = 0$ percepatan konstan sehingga RK4 eksak
  (solusinya polinomial derajat 2), yang dimanfaatkan sebagai uji verifikasi.
- **Lokasi kejadian (event location).** Saat $y$ berganti tanda di dalam satu langkah,
  ukuran sub-langkah $h \in (0, \Delta t]$ dicari dengan bisection sampai presisi mesin
  agar titik tumbukan tepat di $y = 0$. Cara yang sama dipakai untuk titik puncak
  ($v_y = 0$). Hasilnya, jarak, tinggi maksimum, dan waktu tempuh tidak dibatasi
  resolusi $\Delta t$.
- **Presisi ganda.** Semua besaran memakai `real64` dari `iso_fortran_env` (sekitar
  15–16 digit signifikan), jadi galat pembulatan yang terakumulasi selama ribuan
  langkah tetap jauh di bawah galat pemotongan RK4.
- **Waktu dihitung sebagai** $t_n = n\,\Delta t$, bukan dijumlahkan berulang, untuk
  menghindari akumulasi galat pembulatan.

## Struktur direktori

```
fortran-ballistics-sim/
├── src/
│   └── ballistics.f90        # modul ballistics_core + program utama
├── scripts/
│   └── plot_trajectory.py    # visualisasi results/trajectory.csv
├── results/
│   ├── .gitkeep
│   └── trajectory.png        # grafik pratinjau (CSV tidak di-commit)
├── .gitignore
└── README.md
```

## Prasyarat

| Perangkat | Versi | Instalasi |
|---|---|---|
| gfortran | 9 atau lebih baru | Linux: `sudo apt install gfortran` · macOS: `brew install gcc` · Windows: `winget install BrechtSanders.WinLibs.POSIX.UCRT` |
| Python | 3.9 atau lebih baru | `pip install matplotlib pandas numpy` |

## Kompilasi

Dari akar proyek, dengan peringatan ketat dan pemeriksaan runtime:

```bash
gfortran -Wall -Wextra -Wimplicit-interface -fcheck=all -O2 src/ballistics.f90 -o ballistics
```

Untuk produksi (tanpa pemeriksaan runtime, lebih cepat):

```bash
gfortran -O3 -march=native src/ballistics.f90 -o ballistics
```

## Menjalankan

Program harus dijalankan dari akar proyek karena menulis ke `results/trajectory.csv`.

```bash
./ballistics              # Linux/macOS
.\ballistics.exe          # Windows (PowerShell)
```

Contoh keluaran dengan parameter bawaan:

```
==========================================================
  SIMULASI GERAK PELURU 2D  |  RK4, real64
==========================================================
  ...
----------------------------------------------------------
  Jarak terjauh       R     =      105.9576 m
  Tinggi maksimum     H     =       36.8606 m
  Waktu tempuh        T     =        5.4401 s
  Waktu ke puncak     t_H   =        2.4642 s
  Laju tumbukan       v_imp =       25.0551 m/s
  Jumlah titik data         =          5442
----------------------------------------------------------
  Pembanding tanpa hambatan udara (vakum):
  R vakum                   =      254.9291 m
  H vakum                   =       63.7323 m
  Reduksi jarak oleh hambat =       58.4364 %
```

### Format CSV

`results/trajectory.csv` berisi satu baris per langkah waktu dengan 17 digit signifikan.
Baris terakhir adalah titik tumbukan ($y = 0$).

| kolom | arti | satuan |
|---|---|---|
| `t` | waktu | s |
| `x`, `y` | posisi horizontal dan vertikal | m |
| `vx`, `vy` | komponen kecepatan | m/s |

## Visualisasi

```bash
python scripts/plot_trajectory.py
```

Menghasilkan `results/trajectory.png` yang membandingkan lintasan dengan hambatan udara
terhadap parabola vakum dari kondisi awal yang sama. Opsi: `--input`, `--output`,
`--dpi`, dan `--no-vacuum`.

## Verifikasi fisika

```bash
./ballistics --verify
```

Uji yang dijalankan:

1. **$C_d = 0$ vs solusi analitik vakum** untuk $\theta$ = 15°, 30°, 45°, 60°, 75°:

   $$
   R = \frac{v_0^2 \sin 2\theta}{g}, \qquad
   H_{max} = \frac{v_0^2 \sin^2\theta}{2g}, \qquad
   T = \frac{2 v_0 \sin\theta}{g}
   $$

   Lulus bila galat relatif ketiganya ≤ 1×10⁻¹⁰.
2. **Arah efek hambatan.** Dengan $C_d = 0{,}47$, $R$ dan $H$ harus lebih kecil dari nilai vakum.
3. **Konvergensi RK4.** $R$ pada $\Delta t$ = 10⁻² s dan 10⁻³ s harus berselisih relatif ≤ 10⁻⁶.

Hasil pada gfortran 16.2 (Windows, x86-64):

```
  theta     R_sim [m]      R_ref [m]    galat_R    galat_H    galat_T  status
------------------------------------------------------------------------------
   15.0  127.464526622  127.464526622   4.60E-14   1.89E-14   1.63E-14  LULUS
   30.0  220.775036272  220.775036272   3.89E-14   8.75E-14   4.72E-14  LULUS
   45.0  254.929053244  254.929053244   1.86E-13   1.24E-13   8.07E-14  LULUS
   60.0  220.775036272  220.775036272   5.46E-14   1.36E-13   9.53E-14  LULUS
   75.0  127.464526622  127.464526622   2.42E-13   1.41E-13   1.02E-13  LULUS
```

Galat yang tersisa (~10⁻¹³) adalah galat pembulatan `real64`, bukan galat metode.
Program keluar dengan kode 1 bila ada uji yang gagal, sehingga bisa dipakai di CI.

## Mengubah parameter

Semua parameter awal ada di satu tempat, yaitu tipe `sim_params` di bagian atas
`src/ballistics.f90` (blok bertanda `PARAMETER AWAL`). Ubah nilainya lalu kompilasi ulang.

| parameter | bawaan | satuan | keterangan |
|---|---|---|---|
| `v0` | 50.0 | m/s | kecepatan awal |
| `theta_deg` | 45.0 | derajat | sudut elevasi, (−90, 90] |
| `y0` | 0.0 | m | ketinggian awal |
| `mass` | 0.145 | kg | massa peluru |
| `diameter` | 0.074 | m | diameter peluru |
| `cd` | 0.47 | – | koefisien hambat (bola ≈ 0,47) |
| `rho` | 1.225 | kg/m³ | massa jenis udara permukaan laut |
| `g` | 9.80665 | m/s² | gravitasi standar |
| `dt` | 1e-3 | s | langkah waktu RK4 |
| `t_max` | 600 | s | batas waktu pengaman |

Jika `g` diubah, samakan juga konstanta `G` di `scripts/plot_trajectory.py` agar
parabola vakum pembanding tetap benar.

## Penulis

NuRichter
