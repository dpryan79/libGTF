nfail=0
check() {
    out=`diff $1 $2`
    if [ ${#out} -eq 0 ]
    then
        echo pass
    else
        echo fail
        ((nfail++))
    fi
}

#########################################
#Test 1 testFindOverlaps -m any -s ignore
#########################################
echo -n "testFindOverlaps -m any -s ignore ... "
echo -e "A	ENSMUSG00000094121
B	ENSMUSG00000094121
C	ENSMUSG00000094121
D	ENSMUSG00000094121
E	ENSMUSG00000094121
F	ENSMUSG00000094121
G	ENSMUSG00000094121
H	ENSMUSG00000094121
I	ENSMUSG00000094121
J	ENSMUSG00000094121
K	Ambiguous" > correct
./testFindOverlaps -m any -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 2 testFindOverlaps -m exact -s ignore
#########################################
echo -n "testFindOverlaps -m exact -s ignore ... "
echo -e "A	ENSMUSG00000094121
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m exact -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 3 testFindOverlaps -m contain -s ignore
#########################################
echo -n "testFindOverlaps -m contain -s ignore ... "
echo -e "A	ENSMUSG00000094121
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	ENSMUSG00000094121
H	ENSMUSG00000094121
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m contain -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 4 testFindOverlaps -m within -s ignore
#########################################
echo -n "testFindOverlaps -m within -s ignore ... "
echo -e "A	ENSMUSG00000094121
B	ENSMUSG00000094121
C	ENSMUSG00000094121
D	ENSMUSG00000094121
E	ENSMUSG00000094121
F	ENSMUSG00000094121
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	ENSMUSG00000094121
J	ENSMUSG00000094121
K	Ambiguous" > correct
./testFindOverlaps -m within -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 5 testFindOverlaps -m start -s ignore
#########################################
echo -n "testFindOverlaps -m start -s ignore ... "
echo -e "A	ENSMUSG00000094121
B	ENSMUSG00000094121
C	ENSMUSG00000094121
D	ENSMUSG00000094121
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	ENSMUSG00000094121" > correct
./testFindOverlaps -m start -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 6 testFindOverlaps -m end -s ignore
#########################################
echo -n "testFindOverlaps -m end -s ignore ... "
echo -e "A	ENSMUSG00000094121
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	ENSMUSG00000094121
F	ENSMUSG00000094121
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m end -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 7 testFindOverlaps -m any -s same
#########################################
echo -n "testFindOverlaps -m any -s same ... "
echo -e "A	Unassigned_NoFeatures
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	ENSMUSG00000094121
E	Unassigned_NoFeatures
F	ENSMUSG00000094121
G	Unassigned_NoFeatures
H	ENSMUSG00000094121
I	Unassigned_NoFeatures
J	ENSMUSG00000094121
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m any -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 8 testFindOverlaps -m exact -s same
#########################################
echo -n "testFindOverlaps -m exact -s same ... "
echo -e "A	Unassigned_NoFeatures
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m exact -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 9 testFindOverlaps -m contain -s same
#########################################
echo -n "testFindOverlaps -m contain -s same ... "
echo -e "A	Unassigned_NoFeatures
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	ENSMUSG00000094121
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m contain -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 10 testFindOverlaps -m within -s same
#########################################
echo -n "testFindOverlaps -m within -s same ... "
echo -e "A	Unassigned_NoFeatures
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	ENSMUSG00000094121
E	Unassigned_NoFeatures
F	ENSMUSG00000094121
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	ENSMUSG00000094121
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m within -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 11 testFindOverlaps -m start -s same
#########################################
echo -n "testFindOverlaps -m start -s same ... "
echo -e "A	Unassigned_NoFeatures
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	ENSMUSG00000094121
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m start -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 12 testFindOverlaps -m end -s same
#########################################
echo -n "testFindOverlaps -m end -s same ... "
echo -e "A	Unassigned_NoFeatures
B	ENSMUSG00000094121
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	ENSMUSG00000094121
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m end -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 13 testFindOverlaps -m any -s opposite
#########################################
echo -n "testFindOverlaps -m any -s opposite ... "
echo -e "A	ENSMUSG00000094121
B	Unassigned_NoFeatures
C	ENSMUSG00000094121
D	Unassigned_NoFeatures
E	ENSMUSG00000094121
F	Unassigned_NoFeatures
G	ENSMUSG00000094121
H	Unassigned_NoFeatures
I	ENSMUSG00000094121
J	Unassigned_NoFeatures
K	Ambiguous" > correct
./testFindOverlaps -m any -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 14 testFindOverlaps -m exact -s opposite
#########################################
echo -n "testFindOverlaps -m exact -s opposite ... "
echo -e "A	ENSMUSG00000094121
B	Unassigned_NoFeatures
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m exact -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 15 testFindOverlaps -m contain -s opposite
#########################################
echo -n "testFindOverlaps -m contain -s opposite ... "
echo -e "A	ENSMUSG00000094121
B	Unassigned_NoFeatures
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	ENSMUSG00000094121
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m contain -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 16 testFindOverlaps -m within -s opposite
#########################################
echo -n "testFindOverlaps -m within -s opposite ... "
echo -e "A	ENSMUSG00000094121
B	Unassigned_NoFeatures
C	ENSMUSG00000094121
D	Unassigned_NoFeatures
E	ENSMUSG00000094121
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	ENSMUSG00000094121
J	Unassigned_NoFeatures
K	Ambiguous" > correct
./testFindOverlaps -m within -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 17 testFindOverlaps -m start -s opposite
#########################################
echo -n "testFindOverlaps -m start -s opposite ... "
echo -e "A	ENSMUSG00000094121
B	Unassigned_NoFeatures
C	ENSMUSG00000094121
D	Unassigned_NoFeatures
E	Unassigned_NoFeatures
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	ENSMUSG00000094121" > correct
./testFindOverlaps -m start -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 18 testFindOverlaps -m end -s opposite
#########################################
echo -n "testFindOverlaps -m end -s opposite ... "
echo -e "A	ENSMUSG00000094121
B	Unassigned_NoFeatures
C	Unassigned_NoFeatures
D	Unassigned_NoFeatures
E	ENSMUSG00000094121
F	Unassigned_NoFeatures
G	Unassigned_NoFeatures
H	Unassigned_NoFeatures
I	Unassigned_NoFeatures
J	Unassigned_NoFeatures
K	Unassigned_NoFeatures" > correct
./testFindOverlaps -m end -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 19 testCountOverlaps -m any -s ignore
#########################################
echo -n "testCountOverlaps -m any -s ignore ... "
echo -e "A	3
B	3
C	2
D	2
E	3
F	3
G	3
H	3
I	2
J	2
K	4" > correct
./testCountOverlaps -m any -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 20 testCountOverlaps -m exact -s ignore
#########################################
echo -n "testCountOverlaps -m exact -s ignore ... "
echo -e "A	2
B	2
C	0
D	0
E	0
F	0
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m exact -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 21 testCountOverlaps -m contain -s ignore
#########################################
echo -n "testCountOverlaps -m contain -s ignore ... "
echo -e "A	3
B	3
C	0
D	0
E	1
F	1
G	3
H	3
I	0
J	0
K	0" > correct
./testCountOverlaps -m contain -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 22 testCountOverlaps -m within -s ignore
#########################################
echo -n "testCountOverlaps -m within -s ignore ... "
echo -e "A	2
B	2
C	2
D	2
E	2
F	2
G	0
H	0
I	2
J	2
K	4" > correct
./testCountOverlaps -m within -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 23 testCountOverlaps -m start -s ignore
#########################################
echo -n "testCountOverlaps -m start -s ignore ... "
echo -e "A	2
B	2
C	2
D	2
E	0
F	0
G	0
H	0
I	0
J	0
K	2" > correct
./testCountOverlaps -m start -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 24 testCountOverlaps -m end -s ignore
#########################################
echo -n "testCountOverlaps -m end -s ignore ... "
echo -e "A	3
B	3
C	0
D	0
E	3
F	3
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m end -s ignore test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 25 testCountOverlaps -m any -s same
#########################################
echo -n "testCountOverlaps -m any -s same ... "
echo -e "A	0
B	3
C	0
D	2
E	0
F	3
G	0
H	3
I	0
J	2
K	0" > correct
./testCountOverlaps -m any -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 26 testCountOverlaps -m exact -s same
#########################################
echo -n "testCountOverlaps -m exact -s same ... "
echo -e "A	0
B	2
C	0
D	0
E	0
F	0
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m exact -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 27 testCountOverlaps -m contain -s same
#########################################
echo -n "testCountOverlaps -m contain -s same ... "
echo -e "A	0
B	3
C	0
D	0
E	0
F	1
G	0
H	3
I	0
J	0
K	0" > correct
./testCountOverlaps -m contain -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 28 testCountOverlaps -m within -s same
#########################################
echo -n "testCountOverlaps -m within -s same ... "
echo -e "A	0
B	2
C	0
D	2
E	0
F	2
G	0
H	0
I	0
J	2
K	0" > correct
./testCountOverlaps -m within -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 29 testCountOverlaps -m start -s same
#########################################
echo -n "testCountOverlaps -m start -s same ... "
echo -e "A	0
B	2
C	0
D	2
E	0
F	0
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m start -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 30 testCountOverlaps -m end -s same
#########################################
echo -n "testCountOverlaps -m end -s same ... "
echo -e "A	0
B	3
C	0
D	0
E	0
F	3
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m end -s same test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 31 testCountOverlaps -m any -s opposite
#########################################
echo -n "testCountOverlaps -m any -s opposite ... "
echo -e "A	3
B	0
C	2
D	0
E	3
F	0
G	3
H	0
I	2
J	0
K	4" > correct
./testCountOverlaps -m any -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 32 testCountOverlaps -m exact -s opposite
#########################################
echo -n "testCountOverlaps -m exact -s opposite ... "
echo -e "A	2
B	0
C	0
D	0
E	0
F	0
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m exact -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 33 testCountOverlaps -m contain -s opposite
#########################################
echo -n "testCountOverlaps -m contain -s opposite ... "
echo -e "A	3
B	0
C	0
D	0
E	1
F	0
G	3
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m contain -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 34 testCountOverlaps -m within -s opposite
#########################################
echo -n "testCountOverlaps -m within -s opposite ... "
echo -e "A	2
B	0
C	2
D	0
E	2
F	0
G	0
H	0
I	2
J	0
K	4" > correct
./testCountOverlaps -m within -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 35 testCountOverlaps -m start -s opposite
#########################################
echo -n "testCountOverlaps -m start -s opposite ... "
echo -e "A	2
B	0
C	2
D	0
E	0
F	0
G	0
H	0
I	0
J	0
K	2" > correct
./testCountOverlaps -m start -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

#########################################
#Test 36 testCountOverlaps -m end -s opposite
#########################################
echo -n "testCountOverlaps -m end -s opposite ... "
echo -e "A	3
B	0
C	0
D	0
E	3
F	0
G	0
H	0
I	0
J	0
K	0" > correct
./testCountOverlaps -m end -s opposite test.gtf test.sam > seen
check seen correct
rm seen correct

exit $nfail
