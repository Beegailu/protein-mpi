from mpi4py import MPI
import math
from pathlib import Path


# ============================================================
# MEMBACA KOORDINAT ATOM DARI FILE PDB
# ============================================================
def baca_koordinat(nama_file):
    koordinat = []

    with open(nama_file, "r", encoding="utf-8") as f:
        for baris in f:
            if baris.startswith("ATOM"):
                try:
                    x = float(baris[30:38].strip())
                    y = float(baris[38:46].strip())
                    z = float(baris[46:54].strip())

                    koordinat.append((x, y, z))

                except ValueError:
                    continue

    return koordinat


# ============================================================
# MENENTUKAN RANGE PEKERJAAN SETIAP PROCESS
# ============================================================
def tentukan_range(n, rank, size):

    total_pairs = n * (n - 1) // 2

    def pasangan_sebelum(i):
        return i * (2 * n - i - 1) // 2

    target_start = total_pairs * rank // size
    target_end = total_pairs * (rank + 1) // size

    # Mencari indeks awal
    low = 0
    high = n

    while low < high:

        mid = (low + high) // 2

        if pasangan_sebelum(mid) < target_start:
            low = mid + 1
        else:
            high = mid

    start_i = low

    # Mencari indeks akhir
    low = 0
    high = n

    while low < high:

        mid = (low + high) // 2

        if pasangan_sebelum(mid) < target_end:
            low = mid + 1
        else:
            high = mid

    end_i = low

    return start_i, end_i


# ============================================================
# MENGHITUNG INTERAKSI ATOM
# ============================================================
def hitung_interaksi_lokal(
    koordinat,
    start_i,
    end_i,
    threshold
):

    n = len(koordinat)

    jumlah_interaksi = 0
    pasangan_dihitung = 0

    for i in range(start_i, min(end_i, n)):

        xi, yi, zi = koordinat[i]

        for j in range(i + 1, n):

            xj, yj, zj = koordinat[j]

            dx = xi - xj
            dy = yi - yj
            dz = zi - zj

            jarak = math.sqrt(
                dx * dx +
                dy * dy +
                dz * dz
            )

            pasangan_dihitung += 1

            if jarak <= threshold:
                jumlah_interaksi += 1

    return pasangan_dihitung, jumlah_interaksi


# ============================================================
# PROGRAM UTAMA
# ============================================================
def main():

    # --------------------------------------------------------
    # INISIALISASI MPI
    # --------------------------------------------------------
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()


    # --------------------------------------------------------
    # FILE DATASET
    # --------------------------------------------------------
    project_root = Path(__file__).resolve().parents[1]

    nama_file = project_root / "data" / "1LJ4.pdb"

    threshold = 5.0


    # --------------------------------------------------------
    # PROCESS 0 MEMBACA DATA
    # --------------------------------------------------------
    if rank == 0:

        koordinat = baca_koordinat(nama_file)

    else:

        koordinat = None


    # --------------------------------------------------------
    # MEMBAGIKAN DATA KE SEMUA PROCESS
    # --------------------------------------------------------
    koordinat = comm.bcast(
        koordinat,
        root=0
    )

    n = len(koordinat)

    total_pairs = n * (n - 1) // 2


    # --------------------------------------------------------
    # MEMBAGI WORKLOAD
    # --------------------------------------------------------
    start_i, end_i = tentukan_range(
        n,
        rank,
        size
    )


    # Jumlah pasangan yang dikerjakan process ini
    local_pairs = 0

    for i in range(start_i, min(end_i, n)):
        local_pairs += n - i - 1


    # --------------------------------------------------------
    # SINKRONISASI
    # --------------------------------------------------------
    comm.Barrier()

    waktu_mulai = MPI.Wtime()


    # --------------------------------------------------------
    # PERHITUNGAN PARALEL
    # --------------------------------------------------------
    local_pairs_count, local_interactions = (
        hitung_interaksi_lokal(
            koordinat,
            start_i,
            end_i,
            threshold
        )
    )


    # --------------------------------------------------------
    # SELESAI PERHITUNGAN
    # --------------------------------------------------------
    comm.Barrier()

    waktu_selesai = MPI.Wtime()

    local_runtime = waktu_selesai - waktu_mulai


    # --------------------------------------------------------
    # MENGGABUNGKAN HASIL
    # --------------------------------------------------------

    # Total interaksi
    total_interactions = comm.reduce(
        local_interactions,
        op=MPI.SUM,
        root=0
    )


    # Total pasangan
    total_pairs_counted = comm.reduce(
        local_pairs_count,
        op=MPI.SUM,
        root=0
    )


    # Runtime terlama
    runtime_mpi = comm.reduce(
        local_runtime,
        op=MPI.MAX,
        root=0
    )


    # Mengumpulkan jumlah pasangan setiap process
    all_pairs = comm.gather(
        local_pairs,
        root=0
    )


    # --------------------------------------------------------
    # OUTPUT
    # --------------------------------------------------------
    if rank == 0:

        print()
        print("==============================================")
        print("        MPI PROTEIN INTERACTION")
        print("==============================================")

        print(f"Dataset              : data/1LJ4.pdb")
        print(f"Jumlah atom          : {n}")
        print(f"Jumlah pasangan      : {total_pairs}")
        print(f"Threshold Interaksi  : {threshold:.1f} Angstrom")
        print(f"Jumlah process       : {size}")

        print("----------------------------------------------")

        # Menampilkan workload setiap process
        for i in range(size):

            print(
                f"Process {i:<13} : "
                f"{all_pairs[i]} pasangan"
            )

        print("----------------------------------------------")

        print(
            f"Pasangan dihitung    : "
            f"{total_pairs_counted}"
        )

        print(
            f"Jumlah interaksi     : "
            f"{total_interactions}"
        )

        print(
            f"Runtime MPI          : "
            f"{runtime_mpi:.6f} detik"
        )

        print("==============================================")


# ============================================================
# EKSEKUSI
# ============================================================
if __name__ == "__main__":
    main()