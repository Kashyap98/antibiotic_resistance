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
                final_info = data[-1]

                for info in data:
                    if info != final_info:
                        output_file.write(f"{get_string_for_info(info)},")
                    else:
                        output_file.write(f"{get_string_for_info(info)}\n")
