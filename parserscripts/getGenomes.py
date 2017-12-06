import sys

from parsers.genome_extractor import GenomeExtractor

if __name__ == '__main__':
    file_name = ''
    try:
        file_name = sys.argv[1]
    except IOError as e:
        print(e.errno)
        print(e)
    except IndexError:
        file_name = 'data/NC_020879_phage.gb'

    gext = GenomeExtractor(file_name)

    print(gext.find_neighbours(181, 240, 60, 60))
