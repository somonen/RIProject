#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from geopy.distance import distance
import random
import math

class Unit:
    """
    Represents a client/facility in this current problem.
    Stores the name, latitude and longitude of the unit.
    """
    def __init__(self, name: str, latitude: float, longitude: float) -> 'Unit':
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
    def __str__(self) -> str:
        return f"{self.name} {self.latitude} {self.longitude}"

    def __repr__(self) -> str:
        return f"{self.name} {self.latitude} {self.longitude}"

class Facility(Unit):
    """
    Inherits class Unit, enriching it with methods that are used for arbitrarily positioned facilities optimization method.

    Attributes:
        direction: pseudo-random direction in which the facility will move
        best_latitude: latitude at which the solution of the problem had the best value
        best_longitude: longitude at which the solution of the problem had the best value

    Methods:
        move(distance, bounds): move a facility in a direction with a given distance in given bounds
        revert(distance, bounds): move a facility in the opposite direction with a given distance in given bounds
        change_direction(theta = None): if theta is given, the direction at which the facility will move is set to theta,
            otherwise a pseudo-random direction will be picked
        set_best_coord(): saves its own coordinates as the coordinates with the best solution
        teleport(bounds): teleports a facility in a pseudo-random position within the given bounds
        return_to_best_coord(): sets its current position to coordinates with the best solution previously set by set_best_coord()
    """

    
    def __init__(self, name: str, latitude: float, longitude: float) -> 'Facility':
        super().__init__(name, latitude, longitude)
        self.direction = random.uniform(0, 2*math.pi)
        self.best_latitude = latitude
        self.best_longitude = longitude

    def move(self, distance: float, bounds: tuple) -> None:
        self.latitude += distance * math.cos(self.direction)
        self.longitude += distance * math.sin(self.direction)
        self.latitude = np.clip(self.latitude, bounds[0][0], bounds[0][1])
        self.longitude = np.clip(self.longitude, bounds[1][0], bounds[1][1])

    def revert(self, distance: float) -> None:
        self.latitude -= distance * math.cos(self.direction)
        self.longitude -= distance * math.sin(self.direction)
        
    def change_direction(self, theta: float = None) -> None:
        if theta is None:
            self.direction = random.uniform(math.pi/4, 2*math.pi)
        else:
            self.direction = theta

    def set_best_coord(self) -> None:
        self.best_latitude = self.latitude
        self.best_longitude = self.longitude

    def teleport(self, bounds: tuple) -> None:
        self.latitude = random.uniform(bounds[0][0], bounds[0][1])
        self.longitude = random.uniform(bounds[1][0], bounds[1][1])

    def return_to_best_coord(self) -> None:
        self.latitude = self.best_latitude
        self.longitude = self.best_longitude


def generate_units(names: list[str], latitudes: list[str], longitudes: list[str] = None) -> list[Unit]:
    """
    Generates a list of instances of class Unit.
    If longitudes are not passed, it assumes that the coordinates are in the attribute "latitudes" separated by ", "

    Returns:
        list[Unit]
    """

    names = [" ".join(name.split()[1:]) for name in names.splitlines()]
    latitudes = [" ".join(latitude.split()[1:]) for latitude in latitudes.splitlines()]

    if longitudes is None:
        latitudes, longitudes = zip(*[(latitude, longitude) for latitude, longitude in (coords.split(', ') for coords in latitudes)])
    else:
        longitudes = [" ".join(longitude.split()[1:]) for longitude in longitudes.splitlines()]
        
    if not (len(names) == len(latitudes) == len(longitudes)):
        if not (len(latitudes) == len(longitudes)):
            raise ValueError(
                f"Mismatched input data lengths: "
                f"names: {len(names)}, "
                f"latitudes: {len(latitudes)}, "
                f"longitudes: {len(longitudes)}"
            )
        names = names[:len(latitudes)]
        
    return [Unit(name, float(latitude), float(longitude))
            for name, latitude, longitude in zip(names, latitudes, longitudes)]

def extend_bounds(bounds: tuple[tuple[float, float], tuple[float, float]], dist: float = 0):
    """
    Extends bounds by "dist" kilometers in each direction.

    Returns:
        tuple(tuple(float, float), tuple(float, float)): ((southern bound, northern bound), (western bound, eastern bound))
    """
    
    south = distance(kilometers=dist).destination((bounds[0][0], bounds[1][0]), bearing=180)  # move south
    north = distance(kilometers=dist).destination((bounds[0][1], bounds[1][1]), bearing=0)    # move north
    west  = distance(kilometers=dist).destination((bounds[0][0], bounds[1][0]), bearing=270)  # move west
    east  = distance(kilometers=dist).destination((bounds[0][1], bounds[1][1]), bearing=90)   # move east
    
    new_bounds = (
        (south.latitude, north.latitude),
        (west.longitude, east.longitude)
    )
    return new_bounds
    
def calculate_bounds(customers: list[Unit], extend: float = 0) -> tuple[tuple[float, float], tuple[float, float]]:
    """
    Calculates the domain in which the problem will be observed.
    Sets the coordinate of the most northern Unit as the most northern bound, etc. for all other directions.

    Returns:
        tuple(tuple(float, float), tuple(float, float)): ((southern bound, northern bound), (western bound, eastern bound))
    """
    
    bl_corner = (customers[0].latitude, customers[0].longitude)
    tr_corner = (customers[0].latitude, customers[0].longitude)

    for customer in customers:
        if customer.latitude < bl_corner[0]:
            bl_corner = (customer.latitude, bl_corner[1])
        if customer.longitude < bl_corner[1]:
            bl_corner = (bl_corner[0], customer.longitude)
        if customer.latitude > tr_corner[0]:
            tr_corner = (customer.latitude, tr_corner[1])
        if customer.longitude > tr_corner[1]:
            tr_corner = (tr_corner[0], customer.longitude)

    return extend_bounds([(bl_corner[0], tr_corner[0]), (bl_corner[1], tr_corner[1])], extend)

