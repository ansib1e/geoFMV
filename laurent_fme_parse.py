input_dir = r'\\ykr-geo-etl\ETL\FMEServer_SystemShare_Prod\repositories'

import os
import pandas as pd


def find_fmw_files_in_directory(directory):
    """Recursively find all .fmw files in a directory."""
    fmw_files = []

    for foldername, subfolders, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.fmw'):
                full_filepath = os.path.join(foldername, filename)
                fmw_files.append({
                    'File Name': filename,
                    'Full File Path': full_filepath
                })

    return fmw_files


def main(directory, output_excel_path):
    fmw_files = find_fmw_files_in_directory(directory)
    df = pd.DataFrame(fmw_files)

    if not df.empty:
        df.to_excel(output_excel_path, index=False)
        print(f"Results saved to: {output_excel_path}")
    else:
        print("No .fmw files found.")


if __name__ == '__main__':
    input_directory = input_dir
    output_path = "fme_workspaces_ykrgeotel_prod.xlsx"
    main(input_directory, output_path)
