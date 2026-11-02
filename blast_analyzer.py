import csv

def read_blast_results(filename):

    with open(filename, "r", newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file, delimiter="\t")

        if not reader.fieldnames:
            return [], set()

        columns = set(reader.fieldnames)
        results = list(reader)

    return results, columns


def convert_values(hits):
    numeric_fields = ["pident", "qcovs", "evalue", "bitscore"] 

    for hit in hits:
        for field in numeric_fields:
            if not hit.get(field, "").strip():
                raise ValueError(f"Missing value for {field}")

            hit[field] = float(hit[field])

        if not hit.get("sacc", "").strip():
            raise ValueError("Missing accession")

        if not hit.get("stitle", "").strip():
            hit["stitle"] = "No description available"

    return hits


def rank_hits(hits):
    return sorted(
        hits,
        key=lambda hit: (
            hit["evalue"],
            -hit["bitscore"],
            -hit["pident"],
            -hit["qcovs"]
        )
    )


def display_results(hits):

    print("\nBLAST RESULT ANALYZER")
    print("-" * 90)

    print(
        f"{'Rank':<6}"
        f"{'Accession':<12}"
        f"{'Identity':<12}"
        f"{'Coverage':<12}"
        f"{'E-value':<14}"
        f"{'Bit score':<12}"
        f"Description"
    )

    print("-" * 90)

    for rank, hit in enumerate(hits, start=1):
        print(
            f"{rank:<6}"
            f"{hit['sacc']:<12}"
            f"{hit['pident']:<12.1f}"
            f"{hit['qcovs']:<12.1f}"
            f"{hit['evalue']:<14.2g}"
            f"{hit['bitscore']:<12.1f}"
            f"{hit['stitle']}"
        )


def main():
    filename = input("Enter your filename: ").strip()

    try:
        hits, columns = read_blast_results(filename)
    except (FileNotFoundError, IsADirectoryError, PermissionError) as error:
        print(f"Error reading file: {error}")
        return
    except UnicodeDecodeError:
        print("Error: The file could not be read as text.")
        return

    required_columns = {
        "sacc", "pident", "qcovs",
        "evalue", "bitscore", "stitle"
    }

    if not required_columns.issubset(columns):
        missing = required_columns - columns
        print("Missing required columns:", ", ".join(sorted(missing)))
        return

    if not hits:
        print("No BLAST hits found.")
        return

    try:
        hits = convert_values(hits)
    except (ValueError, AttributeError) as error:
        print("Error: Invalid BLAST result data:", error)
        return

    ranked_hits = rank_hits(hits)

    display_results(ranked_hits)

    print("\nBEST HIT")
    print("Accession:", ranked_hits[0]["sacc"])
    print("Description:", ranked_hits[0]["stitle"])


if __name__ == "__main__":
    main()
