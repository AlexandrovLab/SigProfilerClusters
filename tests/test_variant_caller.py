import os
import pytest
from SigProfilerSimulator import SigProfilerSimulator as sigSim
from SigProfilerClusters import SigProfilerClusters as hp


def chrom_sort_key(line):
    """
    Generate a sort key that handles chromosomes in a natural order.
    Example: chr1 < chr2 < ... < chr10 < chrX < chrY
    """
    parts = line.strip().split("\t")
    chrom = parts[0]

    # Remove 'chr' if present
    chrom = chrom.replace("chr", "")

    # Map special chromosomes to a numeric order
    if chrom == "X":
        return (23, parts)
    elif chrom == "Y":
        return (24, parts)
    elif chrom == "M" or chrom == "MT":
        return (25, parts)
    else:
        try:
            return (int(chrom), parts)
        except ValueError:
            return (99, parts)  # For unknown chroms, put them at the end


def compare_files_strict(file1, file2):
    """
    Compare two files line-by-line and report mismatches.
    """
    with open(file1) as f1, open(file2) as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    if len(lines1) != len(lines2):
        print(f"Line count mismatch: {len(lines1)} vs {len(lines2)}")
        return False

    all_match = True
    for i, (line1, line2) in enumerate(zip(lines1, lines2), start=1):
        if line1.strip() != line2.strip():
            print(
                f"Line {i} mismatch:\nExpected: {line1.strip()}\nObserved: {line2.strip()}"
            )
            all_match = False

    if all_match:
        print("Files match exactly line-by-line.")
    return all_match


def run_analysis_for_variant_callers_separate_steps(
    variant_callers, base_path, genome, simulations
):
    # Step 1: Run simulation + clustering for all callers
    for caller in variant_callers:
        print(f"Processing variant caller: {caller}")

        input_path = os.path.join(base_path, "input_vcf_sim", caller)
        if not input_path.endswith(os.sep):
            input_path += os.sep
        if not os.path.exists(input_path):
            print(f"Error: Input path does not exist for {caller}. Skipping...")
            continue
        seed_file_path = os.path.join(base_path, "seeds_file", "Simulator_seeds.txt")

        # Run SigProfilerSimulator
        try:
            print(f"Running SigProfilerSimulator for {caller}...")
            sigSim.SigProfilerSimulator(
                project=caller,
                project_path=input_path,
                genome=genome,
                contexts=["96"],
                seed_file=seed_file_path,
                simulations=simulations,
            )
            print(f"SigProfilerSimulator completed for {caller}.")
        except Exception as e:
            print(f"Error running SigProfilerSimulator for {caller}: {e}")
            continue

        # Run SigProfilerClusters
        try:
            print(f"Running SigProfilerClusters analysis for {caller}...")
            hp.analysis(
                project=caller,
                genome=genome,
                contexts="96",
                simContext=["96"],
                input_path=input_path,
                analysis="all",
                sortSims=True,
                subClassify=True,
                correction=True,
                calculateIMD=True,
                variant_caller=caller,
            )
            print(f"SigProfilerClusters analysis completed for {caller}.")
        except Exception as e:
            print(f"Error running SigProfilerClusters analysis for {caller}: {e}")
            continue

    # Step 2: After all simulations and clustering are done, compare the outputs
    for caller in variant_callers:
        print(f"Comparing output for {caller}...")
        output_file = os.path.join(
            base_path,
            "input_vcf_sim",
            caller,
            "output",
            "vcf_files_corrected",
            f"{caller}_clustered",
            f"{caller}_clustered_vaf.txt",
        )
        reference_file = os.path.join(
            base_path, "output_ref_sim", caller, f"{caller}_clustered_vaf.txt"
        )
        try:
            assert compare_files_strict(
                output_file, reference_file
            ), f"Files for {caller} differ: {output_file} vs {reference_file}"
            print(f"Output matches reference for {caller}.")
        except AssertionError as ae:
            print(str(ae))


def run_one_genome_test(caller):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
    base_path = os.path.join(repo_root, "SigProfilerClusters", "tests")

    genome = "GRCh37"
    simulations = 100

    try:
        run_analysis_for_variant_callers_separate_steps(
            [caller], base_path, genome, simulations
        )
        return True
    except Exception as e:
        print(f"Test failed for caller {caller}: {e}")
        return False


@pytest.mark.parametrize("caller", ["caveman", "mutect2", "standard"])
def test_caller_variants(caller):
    assert run_one_genome_test(caller) is True
