import read
import pickle
import os

def crush_cube(materials, collection_id, cube_dir, cube_format):
    if not os.path.isdir(cube_dir):
        os.mkdir(cube_dir)
    outfile = os.path.join(cube_dir, cube_format.format(collection_id))
    with open(outfile, 'wb') as f:
        pickle.dump(materials, f)


def compress_sam_output_into_a_cube(sam_output_dir, sam_output_format, cube_dir, cube_format, collection_id, year_filter):

    output_files = read.map_sam_output(sam_output_dir, sam_output_format, year_filter)

    materials = read.sam_output(output_files)

    crush_cube(materials, collection_id, cube_dir, cube_format)

    print("Your data has been crushed into a cube\nYou have one hour to remove your cube")

def main():

    sam_output_dir = r"T:\SAM\Outputs\Python"

    sam_output_format = "Eco_(\d+?)_(\d{4})_daily.out"

    cube_dir = r"T:\SAM\Preprocessed\OutputCubes"

    cube_format = r"cube_{}.p"

    collection_id = "mark_twain"

    year_filter = "2010"

    compress_sam_output_into_a_cube(sam_output_dir, sam_output_format, cube_dir, cube_format, collection_id)


if __name__ == "__main__":
    main()
