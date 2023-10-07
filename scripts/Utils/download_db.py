import requests
import os
from tqdm import tqdm


def download_db_file(url, file_path):
    """
    download iphylo_cmd.db from iphylo web server
    url: from which url download
    file_path: to directory
    """
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for 4xx and 5xx status codes

        total_size = int(response.headers.get('content-length', 0))

        with open(file_path, 'wb') as file, tqdm(
                desc=file_path,
                total=total_size,
                unit='B',
                unit_scale=True,
                unit_divisor=1024,
        ) as bar:
            for data in response.iter_content(chunk_size=1024):
                file.write(data)
                bar.update(len(data))

        print("\nDatabase file downloaded successfully.")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download the database file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    #
    # response = requests.get(url)
    # if response.status_code == 200:
    #     # file_path = os.path.join('data', 'iphylo.db')
    #     with open(file_path, 'wb') as file:
    #         file.write(response.content)
    #     print("Database file downloaded successfully.")
    # else:
    #     print("Failed to download the database file. Please check your internet connection.")
    #
    #
    #
    #
