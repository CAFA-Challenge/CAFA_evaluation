from pathlib import Path


if __name__ == "__main__":

    read_directory = "/media/scott/data/cafa3_submissions/ZhangFreddolinoLab/"
    read_author = "Zhang-Freddolino-Lab"
    write_directory = "/media/scott/data/cafa3_submissions/ExamplePredictionsLab1/"
    write_author = "ExamplePredictionsLab1"

    read_directory_obj = Path(read_directory)
    write_directory_obj = Path(write_directory)
    write_directory_obj.mkdir(parents=True, exist_ok=True)

    read_files = read_directory_obj.glob("*.txt")

    for read_file in read_files:
        write_filename = [write_author] + read_file.name.split("_")[1:]
        write_filename = "_".join(write_filename)
        print(f"{read_file} -> {write_filename}")
        write_path = write_directory_obj / write_filename

        with open(read_file, "r") as read_handle, open(write_path, "w") as write_handle:
            write_handle.write(f"AUTHOR: {write_author}\n")
            # skip the first line of the file being read:
            next(read_handle)

            write_handle.writelines(read_handle.readlines())