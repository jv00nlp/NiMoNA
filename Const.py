import numpy as np
import ModelEssentials as ess


# Parameters for variable count of compartments
TOTAL_COLUMNS = ess.get_total_columns_in_population()
COMPARTMENT_COUNT = ess.count_compartments()
COLS_TO_USE = ess.create_tuple_of_col_num(COMPARTMENT_COUNT)
NONE_FUNCTION = "none"

# Constants
NOT_CONTAINS = -1
DEFAULT_DELIMITER = ','

# Paths
PATH_ADJACENCY_CSV = 'CSVs/AdjacencyNRW.csv'
PATH_POPULATIONS_CSV = 'CSVs/PopulationNRW.csv'
PATH_SIR_ADJACENCY = 'CSVs/ModelAdjacency.csv'
PATH_SIR_PLACEHOLDERS = 'CSVs/Placeholder.csv'

PLACEHOLDER_COL = 0
VALUE_COL = 1
CATEGORY_COL = 2

CATEGORY_PLACEHOLDER = 'p'
CATEGORY_VARIABLE = 'v'
CATEGORY_TIME_DEPENDANT = 't'

INDEX_CITY_NAMES = 0
INDEX_POS_X = 1
INDEX_POS_Y = 2

TIME_STEP = 0.1
DAYS = 365

TOTAL_STEPS = int(DAYS / TIME_STEP)
TIME_STEPS = np.linspace(start=0, stop=DAYS, num=TOTAL_STEPS)
ITERATION_STEPS = range(0, TOTAL_STEPS)

ADJACENCY_CITIES_CSV = np.loadtxt(PATH_ADJACENCY_CSV, delimiter=DEFAULT_DELIMITER, skiprows=1,
                                    usecols=ess.create_tuple_of_col_num(77))

# When more than 3 compartments are used, the parameter 'usecols' needs to be changed. City positions always
# need to be in the last two columns in the order x, y
POPULATION_CSV = np.loadtxt(PATH_POPULATIONS_CSV, delimiter=DEFAULT_DELIMITER, usecols=COLS_TO_USE)
POPULATION = POPULATION_CSV[:, 0]
NUMBER_OF_COMPARTMENTS = len(POPULATION_CSV[0])
COMPARTMENTS = range(0, NUMBER_OF_COMPARTMENTS)

ADJACENCY_MATRIX_CITIES = ADJACENCY_CITIES_CSV + np.transpose(ADJACENCY_CITIES_CSV)
NUMBER_OF_CITIES = len(ADJACENCY_MATRIX_CITIES)
CITIES = range(0, NUMBER_OF_CITIES)

# Model subplot parameters
MODEL_NETWORK_ROWS = 7
MODEL_NETWORK_COLUMNS = 11
MODEL_NETWORK_TOTAL = 77
DEFAULT_PLOT_FONTSIZE = 9
SUBPLOT_TUPLE = range(0, MODEL_NETWORK_TOTAL - MODEL_NETWORK_COLUMNS)

# Parameters for network plot
MINIMUM_POPULATION = 0

CONNECTION_COLOR = 'lightgrey'
ANNOTATION_COLOR = 'black'

ANNOTATION_X_OFFSET = -0.01
ANNOTATION_Y_OFFSET = 0.02

Z_ORDER_NETWORK_LINES = 1
Z_ORDER_CITY_POINTS = 2
Z_ORDER_CITY_NAMES = 3
