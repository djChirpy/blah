
duplicates = {}
oldData = ''
id = ''
with open('sortedData.txt', 'r') as sorted:
    for line in sorted:
        fields = line.split()
        id = fields[0].strip() 
        data = fields[1].strip() 
        if id != 'ID':
            if oldData == '':
                oldData = data
                oldID = id           
                continue
            elif data == oldData:
                print(oldID + '\t' + id)
                try:
                    duplicates[oldID].append(id)
                except KeyError:
                    duplicates[oldID] = [id]
                continue
            else:
                oldData = data
                oldID = id
    
        
with open('duplicates.txt', 'w') as outFile:
    for x in duplicates.keys():
        outFile.write(x)
        for entry in duplicates[x]:
            outFile.write('\t' + entry)
        outFile.write('\n')