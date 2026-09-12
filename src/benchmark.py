import subprocess
import sys
import re
from pathlib import Path
import matplotlib.pyplot as plt


# KONFIGURASI
JUMLAH_RUN = 3
JUMLAH_PROCESS = [2, 4, 8]

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SERIAL_PROGRAM = PROJECT_ROOT / "src" / "protein_serial.py"
MPI_PROGRAM = PROJECT_ROOT / "src" / "protein_mpi.py"
RESULTS_DIR = PROJECT_ROOT / "results"

RESULTS_DIR.mkdir(exist_ok=True)



# SERIAL
def jalankan_serial():
    command = [
        sys.executable,
        str(SERIAL_PROGRAM)
    ]

    hasil = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    output = hasil.stdout + hasil.stderr

    # Runtime Serial
    match = re.search(
        r"Runtime Serial\s*:\s*([0-9.]+)\s*detik",
        output
    )
    if not match:
        print("Output serial:")
        print(output)
        raise RuntimeError(
            "Runtime Serial tidak ditemukan."
        )
    return float(match.group(1))



# MPI
def jalankan_mpi(jumlah_process):
    command = [
        "mpiexec",
        "-n",
        str(jumlah_process),
        sys.executable,
        str(MPI_PROGRAM)
    ]

    hasil = subprocess.run(
        command,
        capture_output=True,
        text=True
    )
    output = hasil.stdout + hasil.stderr

    # Runtime MPI 
    match = re.search(
        r"Runtime MPI\s*:\s*([0-9.]+)\s*detik",
        output
    )

    if not match:
        print(f"\nOutput MPI {jumlah_process} process:")
        print(output)

        raise RuntimeError(
            f"Runtime MPI {jumlah_process} process tidak ditemukan."
        )
    return float(match.group(1))



# BENCHMARK
def benchmark():
    hasil = {}
    print("=" * 60)
    print("          BENCHMARK LYSYZOME MPI")
    print("=" * 60)

    # SERIAL
    print("\n[1] SERIAL - 1 PROCESS")
    hasil_serial = []

    for run in range(1, JUMLAH_RUN + 1):
        print(f"Run {run}...", end=" ")

        runtime = jalankan_serial()

        hasil_serial.append(runtime)

        print(f"{runtime:.6f} detik")

    rata_serial = sum(hasil_serial) / len(hasil_serial)

    hasil[1] = rata_serial

    print(f"Average Serial : {rata_serial:.6f} detik")



    # MPI
    for jumlah_process in JUMLAH_PROCESS:
        print(
            f"\n[{jumlah_process}] MPI - "
            f"{jumlah_process} PROCESS"
        )
        hasil_mpi = []

        for run in range(1, JUMLAH_RUN + 1):
            print(f"Run {run}...", end=" ")
            runtime = jalankan_mpi(jumlah_process)
            hasil_mpi.append(runtime)

            print(f"{runtime:.6f} detik")
        rata_mpi = sum(hasil_mpi) / len(hasil_mpi)

        hasil[jumlah_process] = rata_mpi

        print(
            f"Average MPI {jumlah_process} : "
            f"{rata_mpi:.6f} detik"
        )

    return hasil


#HITUNG SPEEDUP & EFFICIENCY
def hitung_performa(hasil):
    runtime_serial = hasil[1]

    data = []

    for process in [1, 2, 4, 8]:
        runtime = hasil[process]
        speedup = runtime_serial / runtime
        efficiency = (speedup / process) * 100

        data.append({
            "process": process,
            "runtime": runtime,
            "speedup": speedup,
            "efficiency": efficiency
        })

    return data

# GRAFIK 1 EXECUTION TIME
def buat_grafik_runtime(data):

    process = [row["process"] for row in data]
    runtime = [row["runtime"] for row in data]

    plt.figure(figsize=(8, 5))

    plt.plot(
        process,
        runtime,
        marker="o",
        linewidth=2
    )

    for x, y in zip(process, runtime):
        plt.annotate(
            f"{y:.4f}s",
            (x, y),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center"
        )

    plt.xlabel("Jumlah Process (P)")
    plt.ylabel("Waktu Komputasi (detik)")
    plt.title("Execution Time vs Jumlah Process")

    plt.xticks(process)
    plt.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()

    file = RESULTS_DIR / "execution_time.png"

    plt.savefig(file, dpi=300)
    plt.close()

    print(f"Grafik runtime disimpan: {file}")



# GRAFIK 2 SPEEDUP
def buat_grafik_speedup(data):

    process = [row["process"] for row in data]
    speedup = [row["speedup"] for row in data]

    #ideal speedup = jumlah processnya
    ideal = process

    plt.figure(figsize=(8, 5))

    plt.plot(
        process,
        speedup,
        marker="s",
        linewidth=2,
        label="MPI Speedup"
    )

    plt.plot(
        process,
        ideal,
        linestyle="--",
        label="Ideal Speedup"
    )

    for x, y in zip(process, speedup):
        plt.annotate(
            f"{y:.2f}x",
            (x, y),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center"
        )

    plt.xlabel("Jumlah Process (P)")
    plt.ylabel("Speedup (T_serial / T_p)")
    plt.title("Speedup vs Jumlah Process")

    plt.xticks(process)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()

    plt.tight_layout()

    file = RESULTS_DIR / "speedup.png"

    plt.savefig(file, dpi=300)
    plt.close()

    print(f"Grafik speedup disimpan: {file}")



# GRAFIK 3bEFFICIENCY
def buat_grafik_efficiency(data):

    process = [row["process"] for row in data]
    efficiency = [row["efficiency"] for row in data]

    plt.figure(figsize=(8, 5))

    plt.plot(
        process,
        efficiency,
        marker="^",
        linewidth=2,
        label="Parallel Efficiency"
    )

    #idealnya efficiency = 100%
    plt.axhline(
        y=100,
        linestyle="--",
        label="100% Ideal Efficiency"
    )

    for x, y in zip(process, efficiency):
        plt.annotate(
            f"{y:.2f}%",
            (x, y),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center"
        )

    plt.xlabel("Jumlah Process (P)")
    plt.ylabel("Efficiency (%)")
    plt.title("Efficiency vs Jumlah Process")

    plt.xticks(process)
    plt.ylim(0, 110)

    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()

    plt.tight_layout()

    file = RESULTS_DIR / "efficiency.png"

    plt.savefig(file, dpi=300)
    plt.close()

    print(f"Grafik efficiency disimpan: {file}")



# HASIL
def tampilkan_hasil(data):

    print("\n")
    print("=" * 65)
    print("                 HASIL BENCHMARK")
    print("=" * 65)

    print(
        f"{'Process':<10}"
        f"{'Runtime':<18}"
        f"{'Speedup':<15}"
        f"{'Efficiency':<15}"
    )

    print("-" * 65)

    for row in data:

        print(
            f"{row['process']:<10}"
            f"{row['runtime']:<18.6f}"
            f"{row['speedup']:<15.2f}x"
            f"{row['efficiency']:<15.2f}%"
        )

    print("=" * 65)


# MAIN
def main():
    hasil = benchmark()

    data = hitung_performa(hasil)

    tampilkan_hasil(data)

    buat_grafik_runtime(data)
    buat_grafik_speedup(data)
    buat_grafik_efficiency(data)

    print("\nBenchmark selesai.")
    print("Semua hasil berada di folder results/")


if __name__ == "__main__":
    main()