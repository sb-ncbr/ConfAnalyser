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

        self.cp = self.cyclopentane
        self.ch = self.cyclohexane
        self.b = self.benzene
        self.o = self.oxane

        self.load_config()
        pass

    def config_file_exist(self):
        return Path(self.get_path()).is_file()

    def load_config(self):
        if not self.config_file_exist():
            self.create_config_template()

        lines = self.read_file()
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

    def get_path(self) -> str:
        return f"./{self.config_file_name}"

    @staticmethod
    def process_line(line: str) -> float:
        return float(line.split(":")[1].replace(",", ".").replace("\n", "").replace("\r", ""))

    def read_file(self) -> list[str]:
        with open(self.get_path(), "r") as file:
            lines = file.readlines()
            return lines

    def create_config_template(self):
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
Tolerance out: 0.3"""

        with open(self.get_path(), "a") as file:
            file.write(template)
