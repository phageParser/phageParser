import time
from parsers.find_accession import PhageFinder


if __name__ == '__main__':
    current_milli_time = lambda: int(round(time.time() * 1000))
    finder = PhageFinder('data/PhagesDB_Data.txt')

    print(finder.findByPhage('Arbiter', 'B2'))
    print(finder.findByPhage('Ariel', 'J'))