import time

from parsers.find_accession import PhageFinder


def current_milli_time():
    return int(round(time.time() * 1000))


if __name__ == '__main__':
    finder = PhageFinder('data/PhagesDB_Data.txt')

    print(finder.find_by_phage('Arbiter', 'B2'))
    print(finder.find_by_phage('Ariel', 'J'))
