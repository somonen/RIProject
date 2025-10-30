#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Util import Unit, Facility, calculate_bounds, extend_bounds
import numpy as np
from haversine import haversine
from itertools import combinations
import random
import math

class AbstractSolver:
    def __init__(self, customers: list[Unit], facilities: list[Unit], p: int) -> 'AbstractSolver':
        self.customers = customers
        self.facilities = facilities
        self.p = p

    def calculate_total_distance(self, current_facilities: list[Unit]) -> float:
    
        total_distance = 0.0
        
        for customer in self.customers:
            total_distance += self.distance_from_nearest_facility(customer, current_facilities)
        return total_distance

    def distance_from_nearest_facility(self, customer: Unit, current_facilities: list[Unit]) -> float:
        
        min_distance = float('inf')
        for facility in current_facilities:
            current_distance = self.calculate_distance(customer, facility)
            min_distance = min(current_distance, min_distance)
        
        return min_distance

    def calculate_distance(self, customer: Unit, facility: Unit) -> float:
        return haversine((customer.latitude, customer.longitude), (facility.latitude, facility.longitude), unit = 'km')

class BruteForceObnPMed(AbstractSolver):

    def solve(self) -> list[list[Unit], float]:    
        
        best_combination = []
        best_total_distance = float('-inf')
        
        for combination in list(combinations(self.facilities, self.p)):
            current_total_distance = self.calculate_total_distance(combination)
            if current_total_distance > best_total_distance:
                best_total_distance = current_total_distance
                best_combination = combination
    
            
        return best_combination, best_total_distance

    def __str__(self) -> str:
        return "Brute force"

    def __repr__(self) -> str:
        return "Brute force"


class VNSObnPMed(AbstractSolver):
    def __init__(self, customers: list[Unit], facilities: list[Unit], p: int, iters: int = 100) -> 'VNSObnPMed':
        super().__init__(customers, facilities, p)
        self.iters = iters

    def solve(self) -> list[list[Unit], float]:
        
        best_solution = random.sample(self.facilities, self.p)
        best_total_distance = self.calculate_total_distance(best_solution)
        
        for _ in range(self.iters):
            for i in range(1, self.p):
                current_facilities = self.variate_facilities(best_solution, i)
                current_total_distance = self.calculate_total_distance(current_facilities)
                if current_total_distance > best_total_distance:
                    best_total_distance = current_total_distance
                    best_solution = current_facilities
        
        
        return best_solution, best_total_distance

    def variate_facilities(self, current_facilities: list[Unit], k: int) -> list[Unit]:
        new_facilities = []
        for i in range(1, k+1):
            new_facility = random.choice(self.facilities)
            while new_facility in current_facilities or new_facility in new_facilities:
                new_facility = random.sample(self.facilities, 1)[0]
            new_facilities.append(new_facility)
        replace_idxs = random.sample(range(len(current_facilities)), k = k)
    
        j = 0
        for i in replace_idxs:
            current_facilities[i] = new_facilities[j]
            j += 1
            
        return current_facilities

    def __str__(self) -> str:
        return "VNS"

    def __repr__(self) -> str:
        return "VNS"

class SAObnPMed(AbstractSolver):
    def __init__(self, customers: list[Unit], facilities: list[Unit], p: int, temp: float = 1, iters: int = 100) -> 'SAObnPMed':
        super().__init__(customers, facilities, p)
        self.temp = temp
        self.iters = iters
        
    def solve(self) -> list[list[Unit], float]:

        best_solution = random.sample(self.facilities, self.p)
        best_total_distance = self.calculate_total_distance(best_solution)
        best_local_solution = best_solution
        best_local_total_distance = best_total_distance
        
        for i in range(1, self.iters):
            current_solution = self.move(best_local_solution)
            current_total_distance = self.calculate_total_distance(current_solution)
    
            if current_total_distance > best_total_distance:
                best_solution = current_solution
                best_total_distance = current_total_distance
                best_local_solution = current_solution
                best_local_total_distance = current_total_distance
            elif random.random() < self.calculate_temperature(i):
                best_local_solution = current_solution
                best_local_total_distance = current_total_distance
    
        return best_solution, best_total_distance

    def calculate_temperature(self, i: int) -> float:
        return self.temp*(np.log(2)/np.log(i + 1))

    def move(self, current_facilities: list[Unit]) -> list[Unit]:

        new_facility = random.choice(self.facilities)
        while new_facility in current_facilities:
            new_facility = random.choice(self.facilities)
    
        replace_idx = random.choice(range(len(current_facilities)))
        current_facilities[replace_idx] = new_facility
    
        return current_facilities

    def __str__(self) -> str:
        return "Simulated annealing"

    def __repr__(self) -> str:
        return "Simulated annealing"

class TSObnPMed(AbstractSolver):
    def __init__(self, customers: list[Unit], facilities: list[Unit], p: int, iters: int = 100, taboo_memory: int = 0) -> 'TSObnPMed':
        super().__init__(customers, facilities, p)
        self.iters = iters
        self.taboo_list = []
        self.taboo_memory = round(len(self.facilities)/2)

    def solve(self) -> list[list[Unit], float]:

        best_solution = random.sample(self.facilities, self.p)
        best_total_distance = self.calculate_total_distance(best_solution)
        best_local_solution = best_solution
        best_local_total_distance = best_total_distance
        self.taboo_list = best_solution
        
        for i in range(1, self.iters):
            current_solution = self.move(best_local_solution)
            current_total_distance = self.calculate_total_distance(current_solution)
    
            if current_total_distance > best_total_distance:
                best_solution = current_solution
                best_total_distance = current_total_distance
                best_local_solution = current_solution
                best_local_total_distance = current_total_distance

            if len(self.taboo_list) > self.taboo_memory:
                self.taboo_list = self.taboo_list[1:]
    
        return best_solution, best_total_distance

    def move(self, current_facilities: list[Unit]) -> list[Unit]:

        new_facility = random.choice(self.facilities)
        while new_facility in current_facilities or new_facility in self.taboo_list:
            new_facility = random.choice(self.facilities)
    
        replace_idx = random.choice(range(len(current_facilities)))
        current_facilities[replace_idx] = new_facility
    
        return current_facilities

    def __str__(self) -> str:
        return "Taboo search"

    def __repr__(self) -> str:
        return "Taboo search"

class APFObnPMed(AbstractSolver):
    def __init__(self, customers: list[Unit], p: int, bounds: tuple, iters: int = 100) -> 'APFObnPMed':
        super().__init__(customers, None, p)
        self.bounds = bounds
        self.iters = iters
        self.distance = math.dist([self.bounds[0][0], self.bounds[1][0]], [self.bounds[0][1], self.bounds[1][1]])/500.0

    def solve(self) -> tuple[list[Unit], float]:

        facilities = [Facility('f', random.uniform(*self.bounds[0]), 
                               random.uniform(*self.bounds[1])) for _ in range(self.p)]
        best_total_distance = 0

        for iteration in range(1, self.iters):
            for facility in facilities:
                current_total_distance = 0
                facility.move(self.distance, self.bounds)
                current_total_distance = self.calculate_total_distance(facilities)
                    
                if current_total_distance > best_total_distance:
                    best_total_distance = current_total_distance
                    facility.set_best_coord()
                else:
                    facility.revert(self.distance)
                    facility.change_direction()

            if random.random() < 0.01:
                for facility in facilities:
                    if random.random() < 1.0/len(facilities):
                        facility.teleport(self.bounds)            
                    
        for facility in facilities:
            facility.return_to_best_coord()
        
        return facilities, best_total_distance

    def __str__(self) -> str:
        return "Arbitrarily positioned"

    def __repr__(self) -> str:
        return "Arbitrarily positioned"


class APFObnPMedWithConstraint(APFObnPMed):
    def __init__(self, customers: list[Unit], p: int, bounds: tuple, min_distance: float, iters: int = 100) -> 'APFObnPMedWithConstraint':
        super().__init__(customers, p, bounds, iters)
        self.bounds = bounds
        self.min_distance = min_distance
        self.iters = iters
        self.distance = math.dist([self.bounds[0][0], self.bounds[1][0]], [self.bounds[0][1], self.bounds[1][1]])/500.0

    def solve(self) -> tuple[list[Unit], float]:

        facilities, best_total_distance = super().solve()

        for facility in facilities:
            facility.return_to_best_coord()
            
            nearest_customer, dist = self.get_nearest_customer(facility)
            i = 0
            while dist < self.min_distance:
                direction = math.atan2(facility.latitude - nearest_customer.latitude, facility.longitude - nearest_customer.longitude)
                if i < 100:
                    facility.change_direction(direction)
                facility.move(self.distance, extend_bounds(self.bounds, self.min_distance))
                nearest_customer, dist = self.get_nearest_customer(facility)
                i += 1
            facility.set_best_coord()
            facility.return_to_best_coord()

        best_total_distance = self.calculate_total_distance(facilities)
        
        return facilities, best_total_distance

    def get_nearest_customer(self, facility: Unit) -> list[Unit, float]:
        nearest_customer = (self.customers[0], self.calculate_distance(self.customers[0], facility))
        for customer in self.customers:
            current_distance = self.calculate_distance(customer, facility)
            if current_distance < nearest_customer[1]:
                nearest_customer = (customer, current_distance)
        
        return nearest_customer

    def __str__(self) -> str:
        return "Arbitrarily positioned"

    def __repr__(self) -> str:
        return "Arbitrarily positioned"

