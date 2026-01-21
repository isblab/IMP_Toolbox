import argparse
import hashlib
import os

def compute_sha256(file_path: str) -> str:
    """Compute the SHA256 checksum of a file.

    Args:
        file_path (str): Path to the file.

    Returns:
        str: SHA256 checksum of the file.
    """
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def verify_sha256(file_path: str, expected_checksum: str) -> bool:
    """Verify the SHA256 checksum of a file.

    Args:
        file_path (str): Path to the file.
        expected_checksum (str): Expected SHA256 checksum.

    Returns:
        bool: True if the checksum matches, False otherwise.
    """
    computed_checksum = compute_sha256(file_path)
    return computed_checksum == expected_checksum

def parse_sha256_file(sha256_file_path: str) -> dict:
    """Parse a SHA256 checksum file.

    Args:
        sha256_file_path (str): Path to the SHA256 checksum file.

    Returns:
        dict: A dictionary mapping file names to their SHA256 checksums.
    """
    checksums = {}
    with open(sha256_file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                checksum, filename = parts
                checksums[filename] = checksum
    return checksums

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify SHA256 checksum of a file.")
    parser.add_argument("--file_path", type=str, help="Path to the file to verify.")
    parser.add_argument("--sha256_file", type=str, default="SHA256SUM", help="Path to the SHA256 checksum file.")

    args = parser.parse_args()

    expected_checksums = parse_sha256_file(args.sha256_file)
    file_name = os.path.basename(args.file_path)

    if verify_sha256(args.file_path, expected_checksums.get(file_name, "")):
        print(f"Checksum verification passed for {args.file_path}.")
        exit(0)
    else:
        print(f"Checksum verification failed for {args.file_path}.")
        exit(1)