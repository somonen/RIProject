#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
from geopy.distance import distance
import random
import math

class Unit:
    def __init__(self, name: str, latitude: float, longitude: float) -> 'Unit':
        self.name = name
        self.coord = np.array([latitude, longitude])
        self.latitude = latitude
        self.longitude = longitude
    def __str__(self) -> str:
        return f"{self.name} {self.coord[0]} {self.coord[1]}"

    def __repr__(self) -> str:
        return f"{self.name} {self.coord[0]} {self.coord[1]}"

class Facility(Unit):
    def __init__(self, name: str, latitude: float, longitude: float) -> 'Facility':
        super().__init__(name, latitude, longitude)
        self.coord = np.array([latitude, longitude])
        self.direction = random.uniform(0, 2*math.pi)
        self.best_coord = self.coord

    def move(self, distance: float, bounds: tuple) -> None:
        self.latitude += distance * math.cos(self.direction)
        self.longitude += distance * math.sin(self.direction)
        self.latitude = np.clip(self.latitude, bounds[0][0], bounds[0][1])
        self.longitude = np.clip(self.longitude, bounds[1][0], bounds[1][1])
        self.coord = np.array([self.latitude, self.longitude])

    def revert(self, distance: float) -> None:
        self.latitude -= distance * math.cos(self.direction)
        self.longitude -= distance * math.sin(self.direction)
        self.coord = np.array([self.latitude, self.longitude])
        
    def change_direction(self, theta: float = None) -> None:
        if theta is None:
            self.direction = random.uniform(math.pi/4, 2*math.pi)
        else:
            self.direction = theta

    def set_best_coord(self) -> None:
        self.best_coord = self.coord

    def teleport(self, bounds: tuple) -> None:
        self.latitude = random.uniform(bounds[0][0], bounds[0][1])
        self.longitude = random.uniform(bounds[1][0], bounds[1][1])

    def return_to_best_coord(self) -> None:
        self.coord = self.best_coord
        self.latitude = self.coord[0]
        self.longitude = self.coord[1]


def generate_units(names: list[str], latitudes: list[str], longitudes: list[str] = None) -> list[Unit]:

    names = [" ".join(name.split()[1:]) for name in names.splitlines()]
    latitudes = [" ".join(latitude.split()[1:]) for latitude in latitudes.splitlines()]

    if longitudes is None:
        latitudes, longitudes = zip(*[(latitude, longitude) for latitude, longitude in (coords.split(', ') for coords in latitudes)])
    else:
        longitudes = [" ".join(longitude.split()[1:]) for longitude in longitudes.splitlines()]
        
    if not (len(names) == len(latitudes) == len(longitudes)):
        raise ValueError(
            f"Mismatched input data lengths: "
            f"names: {len(names)}, "
            f"latitudes: {len(latitudes)}, "
            f"longitudes: {len(longitudes)}"
        )
        
    return [Unit(name, float(latitude), float(longitude))
            for name, latitude, longitude in zip(names, latitudes, longitudes)]

def extend_bounds(bounds: tuple[tuple[float, float], tuple[float, float]], dist: float = 0):
    
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
    bl_corner = (customers[0].coord[0], customers[0].coord[1])
    tr_corner = (customers[0].coord[0], customers[0].coord[1])

    for customer in customers:
        if customer.coord[0] < bl_corner[0]:
            bl_corner = (customer.coord[0], bl_corner[1])
        if customer.coord[1] < bl_corner[1]:
            bl_corner = (bl_corner[0], customer.coord[1])
        if customer.coord[0] > tr_corner[0]:
            tr_corner = (customer.coord[0], tr_corner[1])
        if customer.coord[1] > tr_corner[1]:
            tr_corner = (tr_corner[0], customer.coord[1])

    return extend_bounds([(bl_corner[0], tr_corner[0]), (bl_corner[1], tr_corner[1])], extend)

