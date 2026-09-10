# protein-mpi
TUGAS KTP MPI

1. Deskripsi project

Paralelisasi Perhitungan Interaksi Antaratom pada Struktur Protein 1LJ4 Menggunakan MPI.

2. Struktur project

protein-mpi/
├── src/
│   ├── protein_serial.py
│   └── protein_mpi.py
├── data/
│   └── 1LJ4.pdb
├── results/
└── README.md

3. Persiapan environment

python --version
mpiexec -help

4. Menjalankan Serial

python src/protein_serial.py


5. Menjalankan MPI

mpiexec -n 2 python src/protein_mpi.py
mpiexec -n 4 python src/protein_mpi.py
mpiexec -n 8 python src/protein_mpi.py


6. Hasil pengujian

Masukkan tabel:
Metode	Proses	Runtime
Serial	1	0.549987 s
MPI	2	0.314853 s
MPI	4	0.139034 s
MPI	8	0.116787 s

7. Penjelasan singkat

Serial menjalankan seluruh pasangan secara sequential.
MPI membagi pasangan atom ke beberapa process.
MPI_Bcast digunakan untuk membagikan koordinat.
MPI_Reduce digunakan untuk menggabungkan hasil.
Hasil serial dan MPI sama, yaitu 23.485 interaksi.
