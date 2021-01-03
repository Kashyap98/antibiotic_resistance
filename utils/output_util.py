import os
import utils.dir_utils as dir_utils


def get_string_for_info(info):

    if isinstance(info, list):
        return str(info).replace(",", "+")
    else:
        return str(info)


class OutputFile:

    def __init__(self, file_path=None, header_list=None):
        self.file_path = file_path
        self.header_list = header_list

        self.init_output_file()

    def init_output_file(self):

        final_header = self.header_list[-1]
        with open(self.file_path, "w") as output_file:
            for header in self.header_list:
                if header != final_header:
                    output_file.write(f"{header},")
                else:
                    output_file.write(f"{header}\n")

    def write_data_dict_to_output_file(self, data_dict):
        with open(self.file_path, "a") as output_file:
            for gene, data in data_dict.items():
                output_file.write(f"{gene},")
                final_info_location = len(data)
                count = 0
                for info in data:
                    count += 1
                    if count != final_info_location:
                        output_file.write(f"{get_string_for_info(info)},")
                    else:
                        output_file.write(f"{get_string_for_info(info)}\n")

    def write_data_list_to_output_file(self, data_list):
        final_info_location = len(data_list)
        count = 0
        with open(self.file_path, "a") as output_file:
            for info in data_list:
                count += 1
                if count != final_info_location:
                    output_file.write(f"{get_string_for_info(info)},")
                else:
                    output_file.write(f"{get_string_for_info(info)}\n")
