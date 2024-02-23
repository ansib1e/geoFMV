# geoFMV (Version 1.0)
# Hubert Yoo
# Project Started: September 7, 2023
import math
import os
from datetime import datetime, timedelta

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

# Parses an input video file containing FMV metadata and creates the relavent geojson files needed in STAC format
# to display the FMV overlayed on webmap in the browser
import klvdata
from collections import OrderedDict
import json
import cv2

from geojson import Feature, Polygon, Point, FeatureCollection, dump, dumps, geometry
geometry.DEFAULT_PRECISION = 16

from geojsonio import display
import numpy as np


from json import loads, dumps
from collections import OrderedDict

def to_dict(input_ordered_dict):
    return loads(dumps(input_ordered_dict))


def distribute_frames(total_frames, total_packets):
    # Adjust for the first frame having a metadata packet
    inner_frames = total_frames - 1
    inner_packets = total_packets - 1

    # Calculate the base number of frames for each inner packet
    frames_per_packet = inner_frames // inner_packets

    # Calculate the number of frames that are left after distributing packets
    remaining_frames = inner_frames % inner_packets

    return frames_per_packet, remaining_frames

def parseBinary(input_video, input_bin):
    # Reads input binary file and parses the KLV metadata
    input_filename = os.path.splitext(os.path.basename(input_video))[0]
    intput_filename_wExt = os.path.basename(input_video)


    packet_counter = 0
    stac_bbox = []
    sensor_center_array = []
    frame_center_array = []
    frame_corners_array = []

    cap = cv2.VideoCapture(input_video)
    video_framerate =  cap.get(cv2.CAP_PROP_FPS)
    video_framecount = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    video_framewidth = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    video_frameheight = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))



    with open(input_bin, 'rb') as f:

        for packet in klvdata.StreamParser(f):
            packet_counter += 1
            # extract desired metadata and write to geojson here
            if packet_counter == 1:
                print()
                print("----------------- packet.structure(): ---------------------")
                print(packet.structure())
                # Get BBOX from first packet frame

            #print("---------------- packet.MetdataList(): --------------------")
            packet_metadata = packet.MetadataList()
            print(packet_metadata)
            packet_metadata_json = to_dict(packet_metadata)
            #print("---------------- Packet Metadaata JSON --------------------")
            #print("[Frame:{:3}] ".format(str(frame_counter)) + str(packet_metadata_json))
            """
              Store required metadata (frame corners + center & sensor center)
                'klv_data.misb0601.UserDefinedTimeStamp'
                'klvdata.misb0601.SensorLatitude'
                'klvdata.misb0601.SensorLongitude'
                'klvdata.misb0601.SensorTrueAltitude'
                'klvdata.misb0601.SensorHorizontalFieldOfView'
                'klvdata.misb0601.SensorVerticalFieldOfView'
                'klvdata.misb0601.FrameCenterLatitude'
                'klvdata.misb0601.FrameCenterLongitude'
                'klvdata.misb0601.FrameCenterElevation'
                'klvdata.misb0601.OffsetCornerLatitudePoint1'
                'klvdata.misb0601.OffsetCornerLongitudePoint1'
                'klvdata.misb0601.OffsetCornerLatitudePoint2'
                'klvdata.misb0601.OffsetCornerLongitudePoint2'
                'klvdata.misb0601.OffsetCornerLatitudePoint3'
                'klvdata.misb0601.OffsetCornerLongitudePoint3'
                'klvdata.misb0601.OffsetCornerLatitudePoint4'
                'klvdata.misb0601.OffsetCornerLongitudePoint4'
                
                
            
            currentFrame_TimeStamp = packet_metadata_json['2'][3]
            
            currentFrame_Sensor_Latitude = packet_metadata_json['13'][3]
            currentFrame_Sensor_Longitude = packet_metadata_json['14'][3]
            currentFrame_Sensor_TrueAltitude = packet_metadata_json[2][3]


            currentFrame_Frame_CenterLatitude = packet_metadata_json[2][3]
            currentFrame_Frame_CenterLongitude = packet_metadata_json[2][3]
            currentFrame_Frame_CenterLongitude = packet_metadata_json[2][3]

            currentFrame_Offset_CornerLatitudePoint1 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLongitudePoint1 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLatitudePoint2 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLongitudePoint2 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLatitudePoint3 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLongitudePoint3 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLatitudePoint4 = packet_metadata_json[2][3]
            currentFrame_Offset_CornerLongitudePoint4 = packet_metadata_json[2][3]
            print(currentFrame_TimeStamp)
            """
            currentFrame_TimeStamp = packet_metadata[2][3]
            currentFrame_Sensor_Latitude = float(packet_metadata[13][3])
            currentFrame_Sensor_Longitude = float(packet_metadata[14][3])
            currentFrame_Sensor_TrueAltitude = float(packet_metadata[15][3])
            currentFrame_Sensor_Coords = (currentFrame_Sensor_Longitude, currentFrame_Sensor_Latitude)
            currentFrame_Sensor_Geometry = Point(currentFrame_Sensor_Coords)
            currentFrame_Sensor_Feature = Feature(geometry=currentFrame_Sensor_Geometry,
                                                  properties={"frame_TimeStamp" : currentFrame_TimeStamp,
                                                              "packet_num": packet_counter})


            currentFrame_Frame_CenterLatitude = float(packet_metadata[23][3])
            currentFrame_Frame_CenterLongitude = float(packet_metadata[24][3])
            currentFrame_FrameCenter_Coords = (currentFrame_Frame_CenterLongitude,currentFrame_Frame_CenterLatitude)
            currentFrame_FrameCenter_Geometry = Point(currentFrame_FrameCenter_Coords)
            currentFrame_FrameCenter_Feature = Feature(geometry=currentFrame_FrameCenter_Geometry,
                                                  properties={"frame_TimeStamp" : currentFrame_TimeStamp,
                                                              "packet_num": packet_counter})


            currentFrame_Offset_CornerLatitude_Point1 = float(packet_metadata[26][3])
            currentFrame_Offset_CornerLongitude_Point1 = float(packet_metadata[27][3])
            currentFrame_FrameCorner_TopLeft_Coords = [currentFrame_Frame_CenterLongitude + currentFrame_Offset_CornerLongitude_Point1, currentFrame_Frame_CenterLatitude + currentFrame_Offset_CornerLatitude_Point1]
            currentFrame_Offset_CornerLatitude_Point2 = float(packet_metadata[28][3])
            currentFrame_Offset_CornerLongitude_Point2 = float(packet_metadata[29][3])
            currentFrame_FrameCorner_TopRight_Coords = [currentFrame_Frame_CenterLongitude + currentFrame_Offset_CornerLongitude_Point2, currentFrame_Frame_CenterLatitude + currentFrame_Offset_CornerLatitude_Point2]
            currentFrame_Offset_CornerLatitude_Point3 = float(packet_metadata[30][3])
            currentFrame_Offset_CornerLongitude_Point3 = float(packet_metadata[31][3])
            currentFrame_FrameCorner_BottomLeft_Coords = [currentFrame_Frame_CenterLongitude + currentFrame_Offset_CornerLongitude_Point3, currentFrame_Frame_CenterLatitude + currentFrame_Offset_CornerLatitude_Point3]
            currentFrame_Offset_CornerLatitude_Point4 = float(packet_metadata[32][3])
            currentFrame_Offset_CornerLongitude_Point4 = float(packet_metadata[33][3])
            currentFrame_FrameCorner_BottomRight_Coords = [currentFrame_Frame_CenterLongitude + currentFrame_Offset_CornerLongitude_Point4, currentFrame_Frame_CenterLatitude + currentFrame_Offset_CornerLatitude_Point4]
            currentFrame_FrameCorner_Array = [[currentFrame_FrameCorner_TopLeft_Coords,currentFrame_FrameCorner_TopRight_Coords,currentFrame_FrameCorner_BottomLeft_Coords,currentFrame_FrameCorner_BottomRight_Coords,currentFrame_FrameCorner_TopLeft_Coords]]
            currentFrame_FrameCorner_Geometry = Polygon(currentFrame_FrameCorner_Array)
            currentFrame_FrameCorner_Feature = Feature(geometry=currentFrame_FrameCorner_Geometry,
                                                  properties={"frame_TimeStamp" : currentFrame_TimeStamp,
                                                              "packet_num": packet_counter})
            print("Packet :{:3}".format(packet_counter))
            print(currentFrame_TimeStamp)

            print("------Sensor Center Geometry-------")
            print(currentFrame_Sensor_Feature)
            sensor_center_array.append(currentFrame_Sensor_Feature)

            print("------Frame Center Geometry-------")
            print(currentFrame_FrameCenter_Feature)
            frame_center_array.append(currentFrame_FrameCenter_Feature)

            print("------Frame Corner Geometry-------")
            print(currentFrame_FrameCorner_Feature)
            frame_corners_array.append(currentFrame_FrameCorner_Feature)

            if packet_counter == 1:
                # Get BBOX from first packet frame to for STAC item
                stac_bbox = [currentFrame_Frame_CenterLongitude + currentFrame_Offset_CornerLongitude_Point4, currentFrame_Frame_CenterLatitude + currentFrame_Offset_CornerLatitude_Point4,
                             currentFrame_Frame_CenterLongitude + currentFrame_Offset_CornerLongitude_Point2, currentFrame_Frame_CenterLatitude + currentFrame_Offset_CornerLatitude_Point2]


    # Interpolation
    print("------- Processing: Interpolation ---------")
    # Interpolation Calculations
    interpolation_frames = math.floor(video_framecount/packet_counter) + 1
    print("Interpolated Frames: {}".format(interpolation_frames))
    frame_synced_sensor_center_array = []
    frame_synced_frame_center_array = []
    frame_synced_frame_corners_array = []





    # NO Interpolation, copy over metadata from frames with packets
    for indv_packet in sensor_center_array[:-1]:
        frame_synced_sensor_center_array.extend([indv_packet] * interpolation_frames)

    for indv_packet in frame_center_array[:-1]:
        frame_synced_frame_center_array.extend([indv_packet] * interpolation_frames)

    for indv_packet in frame_corners_array[:-1]:
        frame_synced_frame_corners_array.extend([indv_packet] * interpolation_frames)

    # Interpolate through all packets
    interpolated_sensor_center_array = frame_synced_sensor_center_array
    interpolated_frame_center_array = frame_synced_frame_center_array
    interpolated_frame_corners_array = frame_synced_frame_corners_array

    sensor_center_stack = sensor_center_array
    frame_center_stack = frame_center_array
    frame_corners_stack = frame_corners_array


    # Interpolated to smooth frame transition between packets
    frames_per_packet, remaining_frames = distribute_frames(video_framecount, packet_counter)
    print("[*] Number of additional frames to interpolate between packets: {}".format(str(frames_per_packet)))
    print("[*] Number of frames with no packets remaining at end of video: {}".format(str(remaining_frames)))

    def interpolate_coordinates(start, end, fraction):
        """Interpolate coordinates between two points or polygons."""
        if isinstance(start[0], (list, tuple)):  # Handling polygons
            return [
                interpolate_coordinates(start_coord, end_coord, fraction)
                for start_coord, end_coord in zip(start, end)
            ]
        else:  # Handling points
            return (
                start[0] + fraction * (end[0] - start[0]),
                start[1] + fraction * (end[1] - start[1])
            )

    def interpolate_timestamp(start, end, fraction):
        """Interpolate timestamps."""
        start_time = datetime.fromisoformat(start)
        end_time = datetime.fromisoformat(end)
        delta = end_time - start_time
        interpolated_time = start_time + delta * fraction
        return interpolated_time.isoformat()

    def interpolate_metadata(metadata_array, total_frames, packet_count, geometry_type="Point"):
        """Interpolate metadata for frames."""
        interpolated_features = []

        if geometry_type not in ["Point", "Polygon"]:
            raise ValueError("Supported geometry types are 'Point' and 'Polygon'.")

        # For each known metadata packet, interpolate its position with the next known packet.
        for i in range(packet_count - 1):
            start_feature = metadata_array[i]
            end_feature = metadata_array[i + 1]
            start_coords = start_feature['geometry']['coordinates']
            end_coords = end_feature['geometry']['coordinates']
            start_timestamp = start_feature['properties']['frame_TimeStamp']
            end_timestamp = end_feature['properties']['frame_TimeStamp']

            gap = total_frames // packet_count

            for j in range(gap):
                fraction = j / gap
                coords = interpolate_coordinates(start_coords, end_coords, fraction)
                timestamp = interpolate_timestamp(start_timestamp, end_timestamp, fraction)

                # Create the feature based on the geometry type
                if geometry_type == "Point":
                    geometry = Point(coords)
                else:  # Polygon
                    geometry = Polygon(coords)

                feature = Feature(
                    geometry=geometry,
                    properties={
                        "frame_TimeStamp": timestamp,
                        "packet_num": start_feature['properties']['packet_num'],
                        "frame_number": len(interpolated_features) + 1
                    }
                )
                interpolated_features.append(feature)

        # Append the last known packet data
        last_known = metadata_array[-1].copy()
        last_known['properties']['frame_number'] = len(interpolated_features)
        interpolated_features.append(last_known)

        # If there are any remaining frames, just replicate the last known frame data
        while len(interpolated_features) < total_frames:
            last_feature = interpolated_features[-1].copy()
            last_feature['properties']['frame_number'] = len(interpolated_features)
            interpolated_features.append(last_feature)

        return interpolated_features

    # Test
    #metadata_packets = [...]  # Replace this with the list of your 711 geojson features
    interpolated_sensor_center_array = interpolate_metadata(sensor_center_array, video_framecount, packet_counter, geometry_type="Point")
    interpolated_frame_center_array = interpolate_metadata(frame_center_array, video_framecount, packet_counter, geometry_type="Point")
    interpolated_frame_corners_array = interpolate_metadata(frame_corners_array, video_framecount, packet_counter, geometry_type="Polygon")

    ############################ Creating Feature Collection ###################################

    sensor_center_FeatureCollection = FeatureCollection(frame_synced_sensor_center_array)
    frame_center_FeatureCollection = FeatureCollection(frame_synced_frame_center_array)
    frame_geom_FeatureCollection = FeatureCollection(frame_synced_frame_corners_array)

    # Create Interpolated FeatureCollection
    interpolated_sensor_center_FeatureCollection = FeatureCollection(interpolated_sensor_center_array)
    interpolated_frame_center_FeatureCollection = FeatureCollection(interpolated_frame_center_array)
    interpolated_frame_geom_FeatureCollection = FeatureCollection(interpolated_frame_corners_array)


    print("*"*25)
    # Itterate through FeatureCollections to update properties of each feature/frame to include start and end timestamps
    start_datetime = sensor_center_FeatureCollection['features'][0]["properties"]["frame_TimeStamp"]
    end_datetime = sensor_center_FeatureCollection['features'][-1]["properties"]["frame_TimeStamp"]
    for indv_frame in sensor_center_FeatureCollection["features"]:
        indv_frame["properties"]["start_datetime"] = start_datetime
        indv_frame["properties"]["end_datetime"] = end_datetime
    #print(sensor_center_FeatureCollection)
    #print(frame_center_FeatureCollection)
    #print(frame_geom_FeatureCollection)

    # Save GeoJSON FMV metadata
    output_geojson_sensor_centers_filepath = r'geojson\DayNight\sensor_centers.geojson'
    with open(output_geojson_sensor_centers_filepath, 'w') as f:
        dump(sensor_center_FeatureCollection, f, indent=4)

    output_geojson_frame_centers_filepath = r'geojson\DayNight\frame_centers.geojson'
    with open(output_geojson_frame_centers_filepath, 'w') as f:
        dump(frame_center_FeatureCollection, f, indent=4)

    output_geojson_frame_geoms_filepath = r'geojson\DayNight\frame_geoms.geojson'
    with open(output_geojson_frame_geoms_filepath, 'w') as f:
        dump(frame_geom_FeatureCollection, f, indent=4)

    # Save Interpolated GeoJSON FMV metadata
    output_geojson_sensor_centers_filepath = r'interpolated\DayNight\sensor_centers.geojson'
    with open(output_geojson_sensor_centers_filepath, 'w') as f:
        dump(interpolated_sensor_center_FeatureCollection, f, indent=4)

    output_geojson_frame_centers_filepath = r'interpolated\DayNight\frame_centers.geojson'
    with open(output_geojson_frame_centers_filepath, 'w') as f:
        dump(interpolated_frame_center_FeatureCollection, f, indent=4)

    output_geojson_frame_geoms_filepath = r'interpolated\DayNight\frame_geoms.geojson'
    with open(output_geojson_frame_geoms_filepath, 'w') as f:
        dump(interpolated_frame_geom_FeatureCollection, f, indent=4)
    #print(dumps(frame_geom_FeatureCollection))



    # Create STAC video extension item

    print("-------- Video Properties ----------")
    print("Num. of Packets: {}".format(packet_counter))
    print("Frame Rate: {}".format(str(video_framerate)))
    print("Frame Count: {}".format(str(video_framecount)))
    print("Video Dimensions: ({},{})".format((video_framewidth),str(video_frameheight)))

    # Creating the STAC Item GeoJSON
    STAC_Item_Feature = Feature(properties={
        "title": "STAC_{}".format(input_filename),
        "description": "STAC Item w/ FMV Video Extension Metadata for the '{}' video file.".format(intput_filename_wExt),
        "start_datetime": str(sensor_center_FeatureCollection['features'][0]["properties"]["frame_TimeStamp"]),
        "end_datetime": str(sensor_center_FeatureCollection['features'][-1]["properties"]["frame_TimeStamp"]),
        "proj:epsg": 4326,
        "video:frame_rate": video_framerate,
        "video:frame_count": video_framecount,
        "video:shape": [
          video_framewidth,
          video_frameheight
        ],
        "video:codec_name": "avc1",
        "video:sensor_name": "IR Mitsubishi PtSi Model 500",
        "video:miis_core_id": "0170 F592 F023 7336 4AF8 AA91 62C0 0F2E B2DA 16B7 4341 0008 41A0 BE36 5B5A B96A 3645",
        "datetime": str(sensor_center_FeatureCollection['features'][0]["properties"]["frame_TimeStamp"]),
    },
    geometry=frame_geom_FeatureCollection['features'][0]['geometry'])

    # Add the custom GeoJSON properties specific to the STAC w/ Video Extension
    STAC_Item_Feature["stac_version"] =  "1.0.0"
    STAC_Item_Feature["id"] = "1"
    STAC_Item_Feature["bbox"] = stac_bbox

    STAC_Item_Feature["links"] = [
        {
          "rel": "self",
          "href": "STAC/STAC_Truck.geojson",
          "type": "application/json"
        }
      ]
    STAC_Item_Feature["assets"] = {
        "video": {
          "href": "STAC/" + input_filename + ".mp4",
          "type": "video/mp2t",
          "title": intput_filename_wExt,
          "roles": [
            "data",
            "video",
            "group:video"
          ]
        },
        "video:frame_centers": {
          "title": "frame_centers",
          "href": "STAC/frame_centers.geojson",
          "type": "application/geo+json",
          "roles": [
            "metadata",
            "video:frame_centers",
            "group:video1"
          ]
        },
        "video:frame_geometries": {
          "title": "frame_geometries",
          "href": "STAC/frame_geoms.geojson",
          "type": "application/geo+json",
          "roles": [
            "metadata",
            "video:frame_geometries",
            "group:video1"
          ]
        },
        "video:sensor_centers": {
          "title": "sensor_centers",
          "href": "STAC/sensor_centers.geojson",
          "type": "application/geo+json",
          "roles": [
            "metadata",
            "video:sensor_centers",
            "group:video1"
          ]
        }
      }
    STAC_Item_Feature["stac_extensions"] = [
        "https://stac-extensions.github.io/projection/v1.0.0/schema.json",
        "https://stac-extensions.github.io/video/v1.0.0/schema.json"
      ]
    output_geojson_STAC_filepath = r'geojson\STAC_{}.geojson'.format(input_filename)
    with open(output_geojson_STAC_filepath, 'w') as f:
        dump(STAC_Item_Feature, f, indent=4)
    #print(packet_metadata)
    #print(packet.structure())



def parseFMV():
    print("[+] Completed parsing Full Motion Video file...")
    input_video = r'InputFMV/DayNight/DayFlight.mpg'
    input_bin = r'InputFMV/DayNight/day_flight_mpg_klv_data.out'
    parseBinary(input_video, input_bin)
    return 0

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    parseFMV()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
