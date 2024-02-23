import sys
import klvdata


packet_counter = 0

with open("InputFMV/cheyenne.bin", 'rb') as f:
    print(klvdata.StreamParser(f))
    for packet in klvdata.StreamParser(f):
        print(packet.MetadataList())