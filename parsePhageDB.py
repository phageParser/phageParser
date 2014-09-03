def phageParse(file):
    fulltext = file.read()

    CSV = open('output.CSV', 'w')

    CSV.write('Query, Name, Score, Expect, queryStart, queryEnd, subjectStart, subjectEnd\n')

    #split file up into queries
    queries = fulltext.split('Query=') #queries[0] will be the header
    queries.pop(0)

    for query in queries:
      queryNumber = query.splitlines()[0]
      queryNumber = float(queryNumber)

      #split query into matches
      matches = query.split('>')
      matches.pop(0)

      for match in matches:
        #find name, followed by space or /n
        firstSpace = match.find(' ')
        firstNewline = match.find('\n')
        name = match[ 0 : min(firstSpace, firstNewline) ]

        #find score
        score = match[ match.find('Score = ') + 8 : match.find(' bits ') ]

        #find expect
        expectStart = match.find('Expect = ') + 9
        expect = match[ expectStart : match.find('\n', expectStart) ]

        #find query start
        queryStartStart = match.find('Query: ') + 7
        queryStartEnd = match.find(' ', queryStartStart)
        queryStart = match[ queryStartStart : queryStartEnd ]

        #find query end
        queryEndStart = queryStartEnd
        while match[queryEndStart] == ' ':
          queryEndStart = queryEndStart + 1
        while match[queryEndStart] != ' ':
          queryEndStart = queryEndStart + 1
        queryEndEnd = match.find('\n', queryEndStart)
        queryEnd = int(match[ queryEndStart : queryEndEnd])

        #find sbjct start
        sbjctStartStart = match.find('Sbjct: ') + 7
        sbjctStartEnd = match.find(' ', sbjctStartStart)
        sbjctStart = match[ sbjctStartStart : sbjctStartEnd ]

        #find sbjct end
        sbjctEndStart = sbjctStartEnd
        while match[sbjctEndStart] == ' ':
          sbjctEndStart = sbjctEndStart + 1
        while match[sbjctEndStart] != ' ':
          sbjctEndStart = sbjctEndStart + 1
        sbjctEndEnd = match.find('\n', sbjctEndStart)
        sbjctEnd = int(match[ sbjctEndStart : sbjctEndEnd])

        CSV.write(str(queryNumber) + ', ' + name + ', ' + score + ', ' + expect + ', ' + queryStart + ', ' + str(queryEnd) + ', ' + sbjctStart + ', ' + str(sbjctEnd) + '\n')



file = open('blast-phagesdb.txt', 'r')
phageParse(file)
