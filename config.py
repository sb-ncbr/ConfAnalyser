from pathlib import Path


class CyclohexaneRecord:
    def __init__(self):
        self.t_in = None
        self.t_flat_in = None
        self.t_out = None
        self.t_tw_out = None
        self.t_tw_boat_angle = None
        self.t_angle = None


class CyclopentaneRecord:
    def __init__(self):
        self.t_in = None
        self.t_out = None
        self.t_tw_out = None
        self.t_tw_boat_angle = None
        self.t_angle = None


class BenzeneRecord:
    def __init__(self):
        self.t_flat_in = None


class OxaneRecord:
    def __init__(self):
        self.t_in = None
        self.t_out = None

class OxaneV2Record:
    def __init__(self):
        self.t_flat_angle = None

        self.t_chair_angle = None
        self.t_chair_degree = None

        self.t_boat_angle = None
        self.t_boat_s_degree = None
        self.t_boat_b_degree = None

        self.t_half_angle = None
        self.t_half_s_degree = None
        self.t_half_m_degree = None
        self.t_half_b_degree = None

        self.t_skew_angle = None
        self.t_skew_degree = None

        self.t_envelope_angle = None
        self.t_envelope_s_degree = None
        self.t_envelope_m_degree = None
        self.t_envelope_b_degree = None


class Config:
    """
    Class which loads the properties from the file in the main working directory
    and creates fields to be accessed by classes.
    Contains tolerances for inputs and outputs

    When no config file exists, it creates one with default values
    """

    def __init__(self):
        self.config_file_name = "config.txt"
        self.cyclohexane = CyclohexaneRecord()
        self.cyclopentane = CyclopentaneRecord()
        self.benzene = BenzeneRecord()
        self.oxane = OxaneRecord()
        self.oxanev2 = OxaneV2Record()

        # Aliases to use instead of full names, when one feel brave
        self.cp: CyclopentaneRecord = self.cyclopentane
        self.ch: CyclohexaneRecord = self.cyclohexane
        self.b: BenzeneRecord = self.benzene
        self.o: OxaneRecord = self.oxane
        self.o2: OxaneV2Record = self.oxanev2

        self.load_config()

    def config_file_exist(self):
        """
        Validates whether the file on a set config path exists.
        """
        return Path(self.get_config_file_path()).is_file()

    def load_config(self) -> None:
        """
        Validates if config file exists. IF it does not, create it based on template.
        After that, load the values into their fields.
        """
        if not self.config_file_exist():
            self.create_config_template()

        lines = self.read_file()

        if len(lines) != 41:
            # TODO: Better message + make sure lines count is updated after removal of old oxane code
            print("Config file not loaded properly! Delete the `config.txt` file in the folder with main file and try "
                  "again.")
            exit(-1)

        self.ch.t_in = self.process_line(lines[3])
        self.ch.t_out = self.process_line(lines[4])
        self.ch.t_flat_in = self.process_line(lines[5])
        self.ch.t_tw_out = self.process_line(lines[6])
        self.ch.t_tw_boat_angle = self.process_line(lines[7])
        self.ch.t_angle = self.process_line(lines[8])

        self.cp.t_in = self.process_line(lines[11])
        self.cp.t_out = self.process_line(lines[12])
        self.cp.t_tw_out = self.process_line(lines[13])
        self.cp.t_tw_boat_angle = self.process_line(lines[14])
        self.cp.t_angle = self.process_line(lines[15])

        self.b.t_flat_in = self.process_line(lines[18])

        self.o.t_in = self.process_line(lines[21])
        self.o.t_out = self.process_line(lines[22])

        self.o2.t_flat_angle = self.process_line(lines[25])

        self.o2.t_chair_angle = self.process_line(lines[26])
        self.o2.t_chair_degree = self.process_line(lines[27])

        self.o2.t_boat_angle = self.process_line(lines[28])
        self.o2.t_boat_s_degree = self.process_line(lines[29])
        self.o2.t_boat_b_degree = self.process_line(lines[30])

        self.o2.t_half_angle = self.process_line(lines[31])
        self.o2.t_half_s_degree = self.process_line(lines[32])
        self.o2.t_half_m_degree = self.process_line(lines[33])
        self.o2.t_half_b_degree = self.process_line(lines[34])

        self.o2.t_skew_angle = self.process_line(lines[35])
        self.o2.t_skew_degree = self.process_line(lines[36])

        self.o2.t_envelope_angle = self.process_line(lines[37])
        self.o2.t_envelope_s_degree = self.process_line(lines[38])
        self.o2.t_envelope_m_degree = self.process_line(lines[39])
        self.o2.t_envelope_b_degree = self.process_line(lines[40])

    def get_config_file_path(self) -> str:
        return f"./{self.config_file_name}"

    @staticmethod
    def process_line(line: str) -> float:
        """
        Processes a single line from the config file, taking only float number from
        after the `:` symbol on a given line. Also resolves possible user error
        by using floating point number with `,` instead of `.`
        :return: the loaded floating point number from the line
        """
        return float(line.split(":")[1].replace(",", ".").replace("\n", "").replace("\r", ""))

    def read_file(self) -> list[str]:
        """
        Reads the content of the config file, splits it into 
        list of lines and returns them.
        :return: 
        """
        with open(self.get_config_file_path(), "r") as file:
            lines = file.readlines()
            return lines

    def create_config_template(self) -> None:
        """
        A template config file that will be used when no config is present.
        This function creates the config file with the content of the `template` variable.
        """
        template = """ConfAnalyzer config file. 
Make sure to only edit numeric names and not any name before it.
Cyclohexane:
Tolerance in: 0.1
Tolerance out: 0.6
Tolerance flat in: 0.1
Tolerance twisted out: 0.4
Angle twisted boat: 17.1
Angle tolerance: 1

Cyclopentane:
Tolerance in: 0.1
Tolerance out: 0.6
Tolerance twisted out: 0.54
Angle twisted boat: 10.5
Angle tolerance: 1

Benzene:
Tolerance flat in: 0.1

Oxane:
Tolerance in: 0.1
Tolerance out: 0.3

Oxane V2:
Flat angle tolerance: 5
Chair angle tolerance: 10
Chair expected angle: 35.26
Boat angle tolerance: 10
Boat expected small angle: 35.26
Boat expected big angle: 74.20
Half chair angle tolerance: 10
Half chair expected small angle: 9.07
Half chair expected medium angle: 17.83
Half chair expected big angle: 42.16
Skew angle tolerance: 10
Skew expected angle: 50.84
Envelope angle tolerance: 10
Envelope expected small angle: 17.37
Envelope expected medium angle: 35.26
Envelope expected big angle: 46.86"""

        with open(self.get_config_file_path(), "w") as file:
            file.write(template)
