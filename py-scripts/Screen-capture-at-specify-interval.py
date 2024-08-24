import cv2
import os
from datetime import datetime, timedelta

# PATH:
# C:\Users\User\Desktop\Video-analysis\Pulsatile\Pulsatile.avi
# C:\Users\User\Desktop\Video-analysis\Pulsatile\output
# C:\Users\User\Desktop\Video-analysis\Laminar\Laminar.avi
# C:\Users\User\Desktop\Video-analysis\Laminar\output
# C:\Users\User\Desktop\Video-analysis\Static\Static.avi
# C:\Users\User\Desktop\Video-analysis\Static\output
class VideoFrameCapture:
    def __init__(self, video_file, output_dir='output', capture_interval=1):
        self.video_file = video_file
        self.capture = cv2.VideoCapture(video_file)
        self.output_dir = output_dir
        self.capture_interval = capture_interval  # Interval in seconds

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not self.capture.isOpened():
            raise ValueError(f"Error: Unable to open video file {video_file}")
        else:
            print(f"Successfully opened video file {video_file}")

    def capture_frames(self):
        fps = self.capture.get(cv2.CAP_PROP_FPS)
        if fps == 0:
            raise ValueError(f"Error: Unable to get FPS for video file {self.video_file}")

        frame_count = 0
        next_capture_time = datetime.now()

        while self.capture.isOpened():
            ret, frame = self.capture.read()
            if not ret:
                break  # Break the loop if no more frames are available

            current_time = datetime.now()
            if current_time >= next_capture_time:
                self.save_frame(frame, current_time)
                next_capture_time += timedelta(seconds=self.capture_interval)

            frame_count += 1

        self.capture.release()
        print("Frame capture completed.")

    def save_frame(self, frame, current_time):
        timestamp = current_time.strftime("%Y-%m-%d_%H-%M-%S-%f")[:-3]  # Include milliseconds
        filename = os.path.join(self.output_dir, f"{timestamp}.jpg")
        cv2.imwrite(filename, frame)
        print(f"Saved frame at {filename}")

if __name__ == "__main__":
    video_file_path = input("Enter the path to the video file: ").strip().strip('"')
    output_dir = input("Enter the path to the output directory: ").strip().strip('"')
    capture_interval = float(input("Enter the capture interval in seconds (e.g., 0.25): "))
    
    if not os.path.exists(video_file_path):
        print(f"Error: The video file {video_file_path} does not exist.")
    else:
        try:
            frame_capture = VideoFrameCapture(video_file_path, output_dir, capture_interval)
            frame_capture.capture_frames()
        except ValueError as e:
            print(e)
