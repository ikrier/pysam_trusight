#Function that clips on the left :
def clip_left(read,numcut):
	#Counter for bases to clip in the sequence and quality :
	R=0
	#Counter for bases clipped from the beginning of the alignment :
	M=0
	#process cigar from the left :
	#only read the S from the start : put the "readsoft" as 0 then 1
	readsoft=0
	#Initialize new cigar string :
	newCigar=[]
	#if has "S", count R and remove ; value=4
	for piece in read.cigar:
		if readsoft==0:
			if piece[0]==4:
				R+=piece[1]
	#until M>5:
		#if has "M", count R and remove, count M ; value=0
		#if has "I", count R and remove ; value=1
		#if has "D", count M, remove ; value=2
		#if has "N", count R and remove, count M ; value=3
		#if has "X", count R and remove, count M ; value=8
		if M<numcut:
			if piece[0]==0:
				readsoft=1
				if piece[1]-numcut<=0:
					M+=piece[1]
					R+=piece[1]
				else:
					newCigar.append((0,piece[1]-(numcut-M)))
					R+=(numcut-M)
					M+=(numcut-M)
				continue
			if piece[0]==1:
				readsoft=1
				R+=piece[1]
				continue
			if piece[0]==2:
				readsoft=1
				M+=piece[1]
				continue
			if piece[0]==3:
				readsoft=1
                                if piece[1]-numcut<=0:
                                        M+=piece[1]
                                        R+=piece[1]
                                else:
                                        newCigar.append((3,piece[1]-numcut))
                                        R+=(numcut-M)
					M+=(numcut-M)
				continue
			if piece[0]==8:
				readsoft=1
                                if piece[1]-numcut<=0:
                                        M+=piece[1]
                                        R+=piece[1]
                                else:
                                        newCigar.append((8,piece[1]-numcut))
                                        R+=(numcut-M)
					M+=(numcut-M)
				continue
		if M>=numcut:
			if (len(newCigar)==0) & (piece[0]==2):
				M+=piece[1]
				continue
			else:	
				newCigar.append(piece)
		#process read sequence :
		#remove R bases from the left
	temporaryqual=read.qual[R:]
	read.seq=read.seq[R:]
	#process quality sequence :
		#remove R bases from the left
	read.qual=temporaryqual
	read.cigar=newCigar
	#edit leftmost position
	read.pos+=M

#Function that clips on the right :
def clip_right(read,numcut):
        #Counter for bases to clip in the sequence and quality :
        R=0
        #Counter for bases clipped from the beginning of the alignment :
        M=0
        #process cigar from the left :
        #only read the S from the start : put the "readsoft" as 0 then 1
        readsoft=0
        #Initialize new cigar string :
        newCigar=[]
        #if has "S", count R and remove ; value=4
        for piece in reversed(read.cigar):
                if readsoft==0:
                        if piece[0]==4:
                                R+=piece[1]
        #until M>5:
                #if has "M", count R and remove, count M ; value=0
                #if has "I", count R and remove ; value=1
                #if has "D", count M, remove ; value=2
                #if has "N", count R and remove, count M ; value=3
                #if has "X", count R and remove, count M ; value=8
                if M<numcut:
                        if piece[0]==0:
				readsoft=1
                                if piece[1]-numcut<=0:
                                        M+=piece[1]
                                        R+=piece[1]
                                else:
                                        newCigar.insert(0,(0,piece[1]-(numcut-M)))
                                        R+=(numcut-M)
                                        M+=(numcut-M)
                                continue
                        if piece[0]==1:
				readsoft=1
                                R+=piece[1]
                                continue
                        if piece[0]==2:
				readsoft=1
                                M+=piece[1]
                                continue
                        if piece[0]==3:
				readsoft=1
                                if piece[1]-numcut<=0:
                                        M+=piece[1]
                                        R+=piece[1]
                                else:
                                        newCigar.insert(0,(3,piece[1]-numcut))
                                        R+=(numcut-M)
                                        M+=(numcut-M)
                                continue
                        if piece[0]==8:
				readsoft=1
                                if piece[1]-numcut<=0:
                                        M+=piece[1]
                                        R+=piece[1]
                                else:
                                        newCigar.insert(0,(8,piece[1]-numcut))
                                        R+=(numcut-M)

                                        M+=(numcut-M)
                                continue
                if M>=numcut:
			if (len(newCigar)==0) & (piece[0]==2):
                                M+=piece[1]
                                continue
                        else:
                        	newCigar.insert(0,piece)
                #process read sequence :
                #remove R bases from the left
        temporaryqual=read.qual[:-R]
        read.seq=read.seq[:-R]
        #process quality sequence :
                #remove R bases from the left
        read.qual=temporaryqual
	read.cigar=newCigar
        #don't edit leftmost position


