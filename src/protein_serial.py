import math
import time


def baca_koordinat(data_pdb):
    koordinat = []

    with open(data_pdb, "r") as f:
        for baris in f:
            if baris.startswith("ATOM"):
                x = float(baris[30:38].strip())
                y = float(baris[38:46].strip())
                z = float(baris[46:54].strip())

                koordinat.append((x, y, z))

    return koordinat


def hitung_interaksi_serial(koordinat, threshold):
    n = len(koordinat)

    pasangan_dihitung = 0
    jumlah_interaksi = 0

    #pengukuran waktu dimulai setelah pembacaan file selesai
    waktu_mulai = time.perf_counter()

    for i in range(n):
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

            #innteraksi jika jarak <= threshold
            if jarak <= threshold:
                jumlah_interaksi += 1

    waktu_selesai = time.perf_counter()

    waktu = waktu_selesai - waktu_mulai

    return pasangan_dihitung, jumlah_interaksi, waktu


def main():

    data_pdb = "data/1LJ4.pdb"
    threshold = 5.0

    koordinat = baca_koordinat(data_pdb)

    n = len(koordinat)

    pasangan_dihitung, jumlah_interaksi, waktu = (
        hitung_interaksi_serial(
            koordinat,
            threshold
        )
    )

    print()
    print("==============================================")
    print("        SERIAL PROTEIN INTERACTION")
    print("==============================================")

    print(f"Dataset              : {data_pdb}")
    print(f"Jumlah atom          : {n}")
    print(f"Jumlah pasangan      : {pasangan_dihitung}")
    print(f"Threshold Interaksi  : {threshold:.1f} Angstrom")
    print(f"Jumlah process       : 1")

    print("----------------------------------------------")
    print("HASIL SETIAP PROCESS")
    print("----------------------------------------------")

    print(
        f"Process 0     : "
        f"{pasangan_dihitung:>8} pasangan | "
        f"{jumlah_interaksi:>6} interaksi"
    )

    print("----------------------------------------------")

    print(
        f"Total pasangan       : "
        f"{pasangan_dihitung}"
    )

    print(
        f"Total interaksi      : "
        f"{jumlah_interaksi}"
    )

    print(
        f"Runtime Serial       : "
        f"{waktu:.6f} detik"
    )

    print("==============================================")


if __name__ == "__main__":
    main()