# Dynamic Traffic Assignment (DTA) by [SUMO](https://sumo.dlr.de/index.html)


## Installation

To install and run the calibration tool, follow these steps:

1. Clone the repository to your local machine:

        git clone https://github.com/Khoshkhah/DTA.git     

2. Install [SUMO](https://sumo.dlr.de/docs/Downloads.php).

3. Install the required dependencies. You can use pip to install them:
pip install -r requirements.txt


## Usage

The calibration tool has to be started via [run.py](run.py), which is a command line application. It has the necessary following parameters:

1. -n network file address, that is a sumo network file.

2. -t the trip filename contains trips data. 
        it is a table file contains four columns "vehicle_id", "time", "from_node" and "to_node" that are seperated by comma.

3. -ni the number of iterations.

4. -is the size of each interval in seconds. that is a integer number.

5. -sample-iteration  number of sampling for getting travel time average for each edge, default is 5.

For more information of these input look at the sample grid in the [examples/grid](./grid_sample/) directory.

For getting the other optional arguments use the help command:

         python run.py --help

Also for running the DTA tool, you can use a configuration xml file like [grid.cfg](./grid_sample/grid.cfg):

        python run.py -c grid_sample/grid.cfg


## Contact Information

For any questions, feedback, or inquiries, please feel free to reach out to us:
- Email: kaveh.khoshkhah(at)ut.ee


