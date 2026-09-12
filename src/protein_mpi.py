from mpi4py import MPI
import math
from pathlib import Path

def baca_koordinat(data_pdb):
    koordinat = []
    with open(data_pdb, "r", encoding="utf-8") as f:
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

def tentukan_range(n, rank, size):
    total_pairs = n * (n - 1) // 2

    def pasangan_sebelum(i):
        return i * (2 * n - i - 1) // 2

    target_start = total_pairs * rank // size
    target_end = total_pairs * (rank + 1) // size

    low = 0
    high = n

    while low < high:
        mid = (low + high) // 2

        if pasangan_sebelum(mid) < target_start:
            low = mid + 1
        else:
            high = mid
    start_i = low

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

def hitung_interaksi_lokal(koordinat, start_i, end_i, threshold):
    n = len(koordinat)
    jumlah_interaksi = 0
    pasangan_dihitung = 0

    for i in range(start_i, min(end_i, n)):
        xi, yi, zi = koordinat[i]

        for j in range(i + 1, n):
            xj, yj, zj = koordinat[j]
            #jarak euclidean
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

def main():

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    #lokasi file dataset
    project_root = Path(__file__).resolve().parents[1]
    data_pdb = project_root / "data" / "1LJ4.pdb"
    threshold = 5.0

    #proses 0 membaca dataset
    if rank == 0:
        koordinat = baca_koordinat(data_pdb)
    else:
        koordinat = None
        
    # distribusi data ke semua proses 
    koordinat = comm.bcast(koordinat, root=0)
    n = len(koordinat)
    total_pairs = n * (n - 1) // 2
    
    #pembagian workloadnya 
    start_i, end_i = tentukan_range(
        n,
        rank,
        size
    )

    # hitung jumlah pasangan 
    local_pairs = 0
    for i in range(start_i, min(end_i, n)):
        local_pairs += n - i - 1

    #mulai perhitungan
    comm.Barrier()
    waktu_mulai = MPI.Wtime()
    local_pairs_count, local_interactions = hitung_interaksi_lokal(
        koordinat,
        start_i,
        end_i,
        threshold
    )

    waktu_selesai = MPI.Wtime()
    local_runtime = waktu_selesai - waktu_mulai

    #kumpulan hasil
    total_interactions = comm.reduce(
        local_interactions,
        op=MPI.SUM,
        root=0
    )

    total_pairs_counted = comm.reduce(
        local_pairs_count,
        op=MPI.SUM,
        root=0
    )

    runtime_mpi = comm.reduce(
        local_runtime,
        op=MPI.MAX,
        root=0
    )

    #kumpulan pasangan dan interaksi
    # masing-masing process
    process_results = comm.gather(
        (local_pairs_count, local_interactions),
        root=0
    )

    #output 
    if rank == 0:
        print()
        print("==============================================")
        print("        MPI PROTEIN INTERACTION")
        print("==============================================")

        print(f"Dataset              : data/{data_pdb.name}")
        print(f"Jumlah atom          : {n}")
        print(f"Jumlah pasangan      : {total_pairs}")
        print(f"Threshold Interaksi  : {threshold:.1f} Angstrom")
        print(f"Jumlah process       : {size}")

        print("----------------------------------------------")
        print("HASIL SETIAP PROCESS")
        print("----------------------------------------------")

        for i in range(size):

            pairs, interactions = process_results[i]

            print(
                f"Process {i:<8} : "
                f"{pairs:>8} pasangan | "
                f"{interactions:>6} interaksi"
            )

        print("----------------------------------------------")

        print(
            f"Total pasangan       : "
            f"{total_pairs_counted}"
        )

        print(
            f"Total interaksi      : "
            f"{total_interactions}"
        )

        print(
            f"Runtime MPI          : "
            f"{runtime_mpi:.6f} detik"
        )

        print("==============================================")


if __name__ == "__main__":
    main()