# coding=utf-8
from math import sqrt, floor, ceil   # Mathematische Hilfsfunktionen
from scipy.optimize import fmin      # Zum Ermitteln von Minima einer Funktion
from copy import deepcopy            # Zum Kopieren von Hashtables/Dictionaries
import os                            # Zum Lesen von Ordnern
import datetime                      # Zum Konvertieren von Sekunden in eine Uhrzeit

# Repraesentiert einen Punkt in Form einer Koordinate oder einen Vektor mit dessen x- und y-Wert
class Point:
    def __init__(self, x, y):
        self.x = x      # x-Koordinate
        self.y = y      # y-Koordinate

    # Hash-Funktion
    def __hash__(self):
        return hash((self.x, self.y))

    # Punkt als Zeichenkette
    def __repr__(self):
        return " (" + str(self.x) + " | " + str(self.y) + ") "

    # Vergleich-Funktion
    def __eq__(self, other):
        return (self.x, self.y) == (other.x, other.y)

    # Euklidische Distanz zweier Punkte berechnen (Anwendung des Satz des Pythagoras)
    def euclidean_distance(self, other_point):
        cathetus_a = other_point.x - self.x
        cathetus_b = other_point.y - self.y
        hypotenuse = sqrt(cathetus_a**2 + cathetus_b**2)
        return hypotenuse


# Repraesentiert eine Kante, bzw. eine Strecke mit zwei definierten Punkten A und B
class Edge:
    def __init__(self, point_a, point_b):
        self.point_a = point_a
        self.point_b = point_b

    # Vergleich-Funktion
    def __eq__(self, other):
        return (self.point_a, self.point_b) == (other.point_a, other.point_b) \
               or (self.point_a, self.point_b) == (other.point_b, other.point_a)
    
    # Erhalte Schnittpunkt mit einer anderen Kante, wobei beide Kanten als Geraden fortgefuehrt werden
    def get_intersection_with(self, other_edge):
        # Als allgm. Geradengleichung gilt: y=mx+b.
        # Hier: Zwei Geraden mit jeweils m1 bzwn. m2 und b1 bzw. b2
        # Es gilt: m1 = dy1/dx1 und m2 = dy2/dx2
        dy1 = self.point_b.y - self.point_a.y
        dy2 = other_edge.point_b.y - other_edge.point_a.y
        dx1 = float(self.point_b.x - self.point_a.x)
        dx2 = float(other_edge.point_b.x - other_edge.point_a.x)

        # Wenn Kante a parallel zur y-Achse ist
        if dx1 == 0.0:
            intersection_x = self.point_a.x
            m2 = dy2 / dx2
            b2 = other_edge.point_a.y - m2 * other_edge.point_a.x
            intersection_y = m2 * intersection_x + b2
            return Point(intersection_x, intersection_y)

        # Steigung m1 kann berechnet werden, da dx1 ungleich 0 ist
        m1 = dy1 / dx1
        b1 = self.point_a.y - m1 * self.point_a.x

        # Wenn Kante b parallel zur y-Achse ist
        if dx2 == 0.0:
            intersection_x = other_edge.point_a.x
            intersection_y = m1 * intersection_x + b1
            return Point(intersection_x, intersection_y)

        # Steigung m2 kann berechnet werden, da dx2 ungleich 0 ist
        m2 = dy2 / dx2
        b2 = other_edge.point_a.y - m2 * other_edge.point_a.x
        intersection_x = (b2-b1)/float(m1-m2)
        intersection_y = m1 * intersection_x + b1

        return Point(intersection_x, intersection_y)

    # Prueft, ob Kante sich ausserhalb der Flaeche von zwei aufgespannten Punkten befindet
    def is_out_of_area(self, point_a, point_b):
        # Pruefen, ob Kante im x-Bereich zwischen den Punkten A und B liegt
        if not point_b.x < point_a.x:
            # Punkt A ist links von B oder ueber bzw. unter B
            if self.point_a.x <= point_a.x and self.point_b.x <= point_a.x:
                # Kante liegt links ausserhalb von A
                return True
            if self.point_a.x >= point_b.x and self.point_b.x >= point_b.x:
                # Kante liegt rechts ausserhalb von B
                return True
        if not point_a.x < point_b.x:
            # Punkt A ist rechts von B oder ueber bzw. unter B
            if self.point_a.x >= point_a.x and self.point_b.x >= point_a.x:
                # Kante liegt rechts ausserhalb von A
                return True
            if self.point_a.x <= point_b.x and self.point_b.x <= point_b.x:
                # Kante liegt links ausserhalb von B
                return True

        # Pruefen, ob Kante im y-Bereich zwischen den Punkten A und B liegt
        if not point_b.y < point_a.y:
            # Punkt B ist ueber A oder auf gleicher hoehe
            if self.point_a.y <= point_a.y and self.point_b.y <= point_a.y:
                # Kante liegt unterhalb von A
                return True
            if self.point_b.y >= point_b.y and self.point_a.y >= point_b.y:
                # Kante liegt oberhalb von B
                return True
        if not point_a.y < point_b.y:
            # Punkt A ist ueber B oder auf gleicher hoehe
            if self.point_a.y >= point_a.y and self.point_b.y >= point_a.y:
                # Kante liegt oberhalb von A oder ist auf gleicher hoehe
                return True
            if self.point_a.y <= point_b.y and self.point_b.y <= point_b.y:
                # Kante liegt liegt unterhalb von B oder ist auf gleicher hoehe
                return True

        return False

    # Erhalte den Faktor mit dem die Kante erweitert werden muss, um auf einen bestimmten Punkt zu treffen
    def get_ratio(self, point):
        # Zu pruefende Kante ist parallel zur Y-Achse --> Verhaeltnis ueber y-Werte
        if self.point_b.x == self.point_a.x:
            ratio = (point.y - self.point_a.y) / float(self.point_b.y - self.point_a.y)
        # Sonst: Verhaeltnis ueber x-Werte
        else:
            ratio = (point.x - self.point_a.x) / float(self.point_b.x - self.point_a.x)

        return ratio


# Repraesentiert einen Knoten im Sichtbarkeitsgraphen
class Node:
    def __init__(self, point, id):
        self.point = point  # Jeder Knoten hat einen Punkt
        self.id = id        # Und eine ID, fuer Polygone gilt: ID > 0, Startknoten hat die ID -1, der Endknoten -2

        # Attribute, die fuer den A-Star-Algorithmus relevant sind
        self.g_cost = 0
        self.h_cost = 0
        self.f_cost = 0
        self.previous_node = 0

    # Methode zum Zuweisen der Heuristik (hier: euklidische Distanz)
    def assign_h_cost(self, end_node):
        self.h_cost = self.point.euclidean_distance(end_node.point)
        return self

    # Hash-Funktion, damit Knoten als Schluessel und Werte in einer Hash-Table sein koennen
    def __hash__(self):
        return hash((self.point, self.id))

    # Vergleich-Funktion (Knoten sind identisch, wenn sie den selben Punkt und die selbe ID haben)
    def __eq__(self, other):
        return (self.point, self.id) == (other.point, other.id)

    # Knoten als Zeichenkette
    def __repr__(self):
        if self.id == -2:   # Endknoten
            return " [Y, (" + str(self.point.x) + " | " + str(self.point.y) + ")]"
        if self.id == -1:
            return " [L, (" + str(self.point.x) + " | " + str(self.point.y) + ")]"
        return " [P" + str(self.id) + ", (" + str(self.point.x) + " | " + str(self.point.y) + ")]"

    # Prueft, ob aktueller Knoten A einen anderen Knoten B sehen kann
    def can_see(self, node_b):
        if node_b.id != -2:
            # Zugehoeriges Polygon des ersten Knoten
            polygon_of_first_node = Input.polygon_list[self.id-1]
        # Pruefen ob beide Knoten zu einem Polygon gehoeren
        if self.id == node_b.id:
            for edge in polygon_of_first_node.convert_to_edges():
                # Wenn Knoten sich eine Kante teilen, ...
                current_points = [edge.point_a, edge.point_b]
                if self.point in current_points and node_b.point in current_points:
                    # ... dann sehen sie sich in jedem Fall
                    return True

        # Richtungsvektor der Verbindung zwischen Knoten a und b
        vector_connection = Point(node_b.point.x - self.point.x, node_b.point.y - self.point.y)

        # Potentielle Hindernisse sind jene Kanten, die die Verbindung a zu b lediglich beruehren
        potential_obstacles = []

        # Pruefen, ob eine aller Kanten den Weg von node_a zu node_b verdeckt
        for edge in Input.list_of_all_edges:
            # Richtungsvektor der zu pruefenden Kante
            vector_edge = Point(edge.point_b.x - edge.point_a.x, edge.point_b.y - edge.point_a.y)
            # Pruefen, ob Verbindung a zu b und aktuelle Kante parallel sind
            if vector_connection.x*vector_edge.y == vector_connection.y * vector_edge.x:
                continue    # Kante kann Verbindung nicht schneiden

            # Falls sich die Kante ausserhalb des von a und b aufgespannten Rechtecks befindet
            if edge.is_out_of_area(self.point, node_b.point):
                continue    # Kante kann Verbindung nicht schneiden

            intersection_point = edge.get_intersection_with(Edge(self.point, node_b.point))
            ratio_edge = round(edge.get_ratio(intersection_point), 6)
            ratio_connection = round(Edge(self.point, node_b.point).get_ratio(intersection_point), 6)

            # Kante verdeckt beide Knoten
            if 0 < ratio_edge < 1 and 0 < ratio_connection < 1:
                return False

            # Verbindung koennte durch Eckpunkt eines Polygons gehen, Kanten beruehren die Verbindung a zu b
            elif (ratio_edge == 1 or ratio_edge == 0) and 0 < ratio_connection < 1:
                potential_obstacles.append(edge)

        # Knoten teilt Punkt mit anderen Knoten
        if node_b.point in Input.double_corners or \
                self.point in Input.double_corners:
            return False    # Doppelte Eckpunkte werden ausgelassen

        # Pruefen, ob Verbindung im Polygon liegt bei Knoten gleicher Polygone
        if self.id == node_b.id:
            # OA + AB/2 = OP, OP ist Ortsvektor des zu pruefenden Punktes
            local_vector_checkpoint = Point(self.point.x + vector_connection.x/2,
                                            self.point.y + vector_connection.y/2)
            # Beide Knoten koennen sich nicht sehen, da die Verbindung im Polygon liegt
            if polygon_of_first_node.contains(local_vector_checkpoint):
                return False
            else:
                return True

        # Pruefen, ob Verbindung durch eine Ecke geht, ohne Kanten zu schneiden
        # Dabei muss es benachbarte Kanten geben, die weder Knoten A noch B enthalten
        if len(potential_obstacles) > 1:
            for i in range(len(potential_obstacles)):
                for j in range(i+1, len(potential_obstacles)):
                    points_in_edges = [potential_obstacles[i].point_a, potential_obstacles[i].point_b,
                                       potential_obstacles[j].point_a, potential_obstacles[j].point_b]
                    # Kanten haben einen gemeinsamen Punkt, wenn es Duplikate in der
                    # Liste (bestehend aus jeweils beiden Punkten beider Kanten) gibt
                    if len(points_in_edges) != len(set(points_in_edges)):
                        # Wenn Knoten A und B nicht auf den Kanten liegen
                        if self.point not in points_in_edges and node_b.point not in points_in_edges:
                            return False
        return True
    
    # Prueft, ob Eckpunkt 2 im B-Bereich von Eckpunkt 1 liegt, wodurch die Verbindung von E2 zu E1 redundant waere
    def is_redundant_to(self, node_2):

        # Polygon des Eckpunkts 1
        polygon = Input.polygon_list[self.id-1]

        # Benachbarte Ecken von Eckpunkt 1
        adjacent_points = polygon.get_adjacent_points(self.point)

        # Richtungsvektoren der beiden fortgefuehrten Kanten an Eckpunkt 1
        r1 = Point(2*self.point.x-adjacent_points[0].x, 2*self.point.y-adjacent_points[0].y)
        r2 = Point(2*self.point.x-adjacent_points[1].x, 2*self.point.y-adjacent_points[1].y)

        # Richtungsvektor vom Eckpunkt 1 zum Eckpunkt 2
        vector_connection = Point(node_2.point.x - self.point.x, node_2.point.y - self.point.y)

        # Determinante der Orientierungs-Matrix berechnen
        def det(p, q, r):
            return (q.x * r.y + p.x * q.y + p.y * r.x) - (p.y * q.x + q.y * r.x + p.x * r.y)

        determinant_r1 = det(self.point, Point(r1.x + vector_connection.x, r1.y + vector_connection.y), node_2.point)
        determinant_r2 = det(self.point, Point(r2.x + vector_connection.x, r2.y + vector_connection.y), node_2.point)

        # Falls Eckpunkt 2 zwischen den beiden fortgefuehrten Kanten von Eckpunkt 1 liegt
        # Wenn also die Determinante der Orientierungsmatrix einmal negativ und einmal positiv ist
        if determinant_r1 < 0 < determinant_r2 or determinant_r2 < 0 < determinant_r1:
            return True

        return False


# Repraesentiert ein Polygon mit seiner ID und seinen Eckpunkten
class Polygon:
    def __init__(self, points, id):
        self.points = points
        self.id = id

    # Prueft, ob ein Punkt im Polygon befindet mithilfe der Strahl-Methode
    def contains(self, point):
        polygon_edges = self.convert_to_edges()  # Alle Kanten des Polygons
        count_to_left_of_point = 0               # Variable zum Zaehlen der Schnittpunkte des Strahls links vom Punkt
        for edge in polygon_edges:
            if edge.point_a.y <= point.y <= edge.point_b.y or edge.point_b.y <= point.y <= edge.point_a.y:
                # Punkt ist auf einer Y-Ebene mit der aktuellen Kante des Polygons
                if edge.point_a.y == edge.point_b.y == point.y:
                    # Kante liegt auf Strahl --> Wird aber als darueber gezaehlt
                    continue
                intersection_point = edge.get_intersection_with(Edge(Point(0, point.y), Point(1, point.y)))
                if intersection_point in [edge.point_a, edge.point_b]:
                    # Strahl geht durch Eckpunkt der Kante, Punkte auf dem Strahl werden als Punkte darueber gezaehlt
                    if edge.point_a.y == point.y:
                        # Punkt A der Kante liegt auf dem Strahl
                        if edge.point_b.y > point.y:
                            # Kante schneidet nicht Strahl, da Punkt B und Punkt A "ueber" dem Strahl
                            continue
                    else:
                        # Punkt B liegt auf dem Strahl
                        if edge.point_a.y > point.y:
                            # Kante schneidet nicht Strahl, da Punkt B und Punkt A "ueber" dem Strahl
                            continue

                if intersection_point.x < point.x:
                    # Schnittpunkt der Kante links vom zu pruefenden Punkt
                    count_to_left_of_point += 1

        # Wenn Anzahl der Schnittpunkte links vom Punkt ungerade ist --> Punkt liegt im Polygon
        if count_to_left_of_point % 2 != 0:
            return True
        else:
            return False

    # Erhalte benachbarten Punkte eines Punktes in einem Polygon
    def get_adjacent_points(self, point):
        adjacent_points = []
        for edge in self.convert_to_edges():
            if point == edge.point_a:
                adjacent_points.append(edge.point_b)
            if point == edge.point_b:
                adjacent_points.append(edge.point_a)

        return adjacent_points

    # Zusammenfuegen von Punkten eines Polygons zu Kanten
    def convert_to_edges(self):
        edges_list = []
        number_of_corners = len(self.points)
        for i in range(number_of_corners):
            if i == number_of_corners - 1:
                edges_list.append(Edge(self.points[i], self.points[0]))
            else:
                edges_list.append(Edge(self.points[i], self.points[i + 1]))
        return edges_list

    # Gibt nur konvexe Punkte eines Polygons zurueck
    def convex_vertices_only(self):
        points = self.points
        # Finde oberste linke Ecke, da diese in jedem Fall konvex ist
        top_left_corner_index = 0
        for i in range(1, len(points)):
            if points[i].x <= points[top_left_corner_index].x:
                if points[i].y >= points[top_left_corner_index].y:
                    top_left_corner_index = i

        # neue Liste der Eckpunkte mit der obersten linken Ecke als Startwert
        new_points_list = points[top_left_corner_index:]
        new_points_list.extend(points[0:top_left_corner_index])
        points = new_points_list

        # Alle Ecken, bei denen die Kante rechts von der vorherigen Kante angelegt wird, sind in list_right
        list_right = []
        # Die restlichen Ecken kommen in list_left
        list_left = []

        # Determinante der Orientierungs-Matrix berechnen
        def det(p, q, r):
            return (q.x * r.y + p.x * q.y + p.y * r.x) - (p.y * q.x + q.y * r.x + p.x * r.y)

        # Richtungen fuer jede Ecke bestimmen mithilfe von drei Punkten, der jeweils anliegenden zwei Kanten
        for i in range(len(points)):
            if i == 0:
                point_p = points[len(points)-1]
            else:
                point_p = points[i-1]

            point_q = points[i]

            if i == len(points)-1:
                point_r = points[0]
            else:
                point_r = points[i+1]

            determinant = det(point_p, point_q, point_r)

            # Orientierung ist im Uhrzeigersinn --> Rechtsrum
            if determinant < 0:
                list_right.append(point_q)
            # Orientierung gegen den Uhrzeigersinn --> Linksrum
            else:
                list_left.append(point_q)

        if points[0] in list_right:
            # Ecken, deren Kanten einen Rechts-Turn machen sind konvex

            return list_right
        else:
            # Ecken, deren Kanten einen Links-Turn machen sind konvex
            return list_left


class VisibilityGraph:

    # Gibt einen Sichtbarkeitsgraph als Hash-Table zurueck (Rekursive Methode)
    @staticmethod
    def get_graph(to_visit, current_graph):

        # Es muessen keine Knoten mehr besucht werden
        if len(to_visit) == 1:
            return current_graph

        # Knoten wird abgearbeitet
        current_node = to_visit[0]
        to_visit.pop(0)

        # Jeden zu besuchenden Knoten mit aktuellem Knoten auf Sichtbarkeit und Redundanz pruefen
        for node in to_visit:
            if node.can_see(current_node):
                if not node.is_redundant_to(current_node):
                    current_graph[current_node].append(node)
                if not current_node.is_redundant_to(node):
                    current_graph[node].append(current_node)

        # Naechsten zu besuchenden Knoten pruefen, Erweiterung des aktuellen Sichtbarkeitsgraphen
        return VisibilityGraph.get_graph(to_visit, current_graph)

    # Methode zum Hinzufuegen eines Start- oder Endknotens in einen bereits vorhandenen Sichtbarkeitsgraphen
    @staticmethod
    def add_node(to_add, visibility_graph):
        all_key_nodes = visibility_graph.keys()
        visibility_graph[to_add] = []

        # Jeden Knoten mit zu ergaenzendem Knoten auf Sichtbarkeit und Redundanz pruefen
        for node in all_key_nodes:
            if node.can_see(to_add):
                    if node.id >= 0:
                        if not node.is_redundant_to(to_add):
                            visibility_graph[node].append(to_add)
                            visibility_graph[to_add].append(node)
                    else:
                        visibility_graph[node].append(to_add)
                        visibility_graph[to_add].append(node)

        return visibility_graph


# Klasse, die alle wichtigen Methoden zum A* Pathfinding-Algorithmus enthaelt
class AStar:

    # Findet den kuerzesten Pfad in einem Sichtbarkeitsgraphen
    @staticmethod
    def get_shortest_path(visibility_graph):
        start_node = Input.start_node
        open_list = [start_node]
        current_node = start_node
        closed_list = []

        # Solange die Open-Liste nicht leer ist
        while len(open_list) != 0:

            # Alle sichtbaren Knoten des aktuellen Knotens
            visible_nodes = visibility_graph[current_node]

            # Knoten wird als abgeschlossen markiert
            closed_list.append(current_node)
            open_list.remove(current_node)

            # Iteriere durch sichtbare Knoten des aktuellen Knotens
            for i in range(len(visible_nodes)):
                visible_node = visible_nodes[i]
                # Berechne gCost (Abstand zum Startknoten)
                new_g_cost = current_node.g_cost + current_node.point.euclidean_distance(visible_node.point)

                # Wenn sichtbarer Knoten als abgeschlossen markiert und der neue Pfad zu ihm nicht guenstiger ist
                if visible_node in closed_list and new_g_cost >= visibility_graph[current_node][i].g_cost:
                    continue

                # Falls Knoten noch nicht entdeckt oder ein guenstigerer Pfad gefunden wurde
                if visible_node not in open_list or new_g_cost < visibility_graph[current_node][i].g_cost:
                    # Aktualisieren der gCost und Gesamtkosten fuer diesen Knoten
                    visibility_graph[current_node][i].g_cost = new_g_cost
                    visibility_graph[current_node][i].f_cost = new_g_cost + visible_node.h_cost

                    # Setzen des Zeigers auf den vorherigen Knoten
                    visibility_graph[current_node][i].previous_node = current_node

                    # Falls Knoten noch nicht entdeckt
                    if visible_node not in open_list:
                        # Knoten wird als offen markiert
                        open_list.append(visible_node)

            # Erhalte Knoten mit geringsten Gesamtkosten
            if len(open_list) > 0:
                best_next_node = open_list[0]
                for node in open_list:
                    if node.f_cost < best_next_node.f_cost:
                        best_next_node = node
                current_node = best_next_node

            else:
                return None     # Kein Pfad konnte gefunden werden

            # Wenn End-Knoten Knoten mit geringsten Gesamtkosten ist
            if current_node.id == -2:
                # Optimaler Pfad wurde gefunden, muss rekonstruiert werden
                return AStar.reconstruct_path(current_node)

        # Es konnte kein Pfad gefunden werden
        return None

    # Methode zur Rekonstruierung des optimalsten Wegs
    @staticmethod
    def reconstruct_path(end_node):
        shortest_path = [end_node]
        previous_node = end_node.previous_node
        total_length = end_node.point.euclidean_distance(previous_node.point)

        # Pfad erweitern mithilfe des Vorgaenger-Knotens, Erhoehung der Gesamtlaenge
        while previous_node.id != -1:                   # Solange der Startknoten nicht erreicht ist
            shortest_path.append(previous_node)
            total_length += previous_node.point.euclidean_distance(previous_node.previous_node.point)
            if previous_node.previous_node.id == -1:    # Falls der vorherige Knoten der Startknoten ist
                total_length += previous_node.point.euclidean_distance(previous_node.previous_node.point)
            previous_node = previous_node.previous_node

        shortest_path.append(previous_node)
        return shortest_path, total_length    # Rueckgabe des kuerzesten Wegs und der Gesamtlaenge als Tupel


# Klasse, die allgemeine Methoden zur Loesung der Aufgabe enthaelt
class LisaRennt:

    # Gibt den optimalsten Pfad mit Sichtbarkeitsgraphen zurueck,
    # bei dem Lisa ihr Haus so spaet wie moeglich verlassen muss
    @staticmethod
    def get_optimal_result():
        to_visit = []       # Liste aller Knoten, die auf Sichtbarkeit ueberprueft werden muss
        start_graph = {}    # Anfangsgraph mit allen Polygonecken als Knoten

        # Konvexe Ecken aller Polygone als Ausgangsknoten hinzufuegen
        for polygon in Input.polygon_list:
            for convex_vertice in polygon.convex_vertices_only():
                # Erstellung eines Knotens aus einem konvexen Eckpunkt
                converted_node = Node(convex_vertice, polygon.id)

                # Heuristik fuer Knoten berechnen und hinzufuegen
                converted_node = converted_node.assign_h_cost(Input.end_node)

                # Knoten jeweils initialisieren
                to_visit.append(converted_node)
                start_graph[converted_node] = []

        # Liste aller Kanten generieren (Zur Erstellung des Sichtbarkeitsgraphen erforderlich)
        for polygon in Input.polygon_list:
            Input.list_of_all_edges.extend(polygon.convert_to_edges())

        # Sichtbarkeitsgraph mit Start- aber ohne Endknoten
        visibility_graph = VisibilityGraph.get_graph(to_visit, start_graph)
        visibility_graph = VisibilityGraph.add_node(Input.start_node, visibility_graph)

        # Untersten und obersten Knoten ermitteln (um Bereich zu definieren), Startknoten als Referenzwert
        smallest_y = Input.start_node.point.y
        biggest_y = Input.start_node.point.y
        for node in visibility_graph.keys():
            if node.point.y < smallest_y:
                smallest_y = node.point.y
            if node.point.y > biggest_y:
                biggest_y = node.point.y

        # Maximum von t(x), s. Dokumentation
        optimum_y = fmin(lambda x: -LisaRennt.t(x), 0, disp=False)

        # Bereich justieren, indem optimaler Weg gesucht wird, da das Optimum ein Teil dieses Bereichs sein muss
        if optimum_y < smallest_y:
            smallest_y = int(floor(optimum_y))
        if optimum_y > biggest_y:
            biggest_y = int(ceil(optimum_y))

        # Ermitteln des optimalsten Wegs durch Ausprobieren aller moeglichen Endknoten im vordefinierten Bereich
        full_visibility_graph = VisibilityGraph.add_node(Node(Point(0,smallest_y), -2), deepcopy(visibility_graph))
        best_visibility_graph = deepcopy(full_visibility_graph)
        calc_path = AStar.get_shortest_path(full_visibility_graph)
        if calc_path is None:
            return  # Kein Pfad
        best_path = calc_path
        Output.highest_time_after_departure = LisaRennt.time(smallest_y, calc_path[1])
        for i in range(smallest_y+1, biggest_y):
            full_visibility_graph = VisibilityGraph.add_node(Node(Point(0, i), -2), deepcopy(visibility_graph))
            calc_path = AStar.get_shortest_path(full_visibility_graph)
            time_after_departure = LisaRennt.time(i, calc_path[1])

            if time_after_departure > Output.highest_time_after_departure:
                Output.highest_time_after_departure = time_after_departure

                best_path = calc_path
                best_visibility_graph = deepcopy(full_visibility_graph)

        return best_path, best_visibility_graph  # Rueckgabe des optimalsten Wegs (inkl. Laenge) und Sichtbarkeitsgraph

    # Zeitdifferenz ohne Hindernisse
    @staticmethod
    def t(x):
        # Optimaler Punkt (0|x) auf der y-Achse, wenn es keine Hindernisse gibt
        s_x = Input.start_node.point.x  # Standort Lisa, x-Wert
        s_y = Input.start_node.point.y  # Standort Lisa, y-Wert
        return (3*x/25.0) - (6 * sqrt(s_x**2 + (s_y-x)**2))/25.0

    # Zeitdifferenz mit Hindernissen
    @staticmethod
    def time(x, path_length):
        return (3*x/25.0) - (6*path_length)/25.0


class Input:
    polygon_list = []
    list_of_all_edges = []
    start_node = Node(Point(0,0),-1)     # Platzhalter Start-Knoten
    end_node = Node(Point(0,0), -2)      # Platzhalter End-Knoten
    double_corners = []

    @staticmethod
    def polygon_liste_einlesen(file_name):
        lines = [line.rstrip('\n') for line in open('beispieldaten/' + file_name)]
        count_polygons = int(lines[0])
        all_vertices = []
        for i in range(count_polygons):
            corner_info = lines[i + 1].split(" ")
            count_corners = int(corner_info[0])
            corner_info.pop(0)
            corner_list = []
            x = 0
            for j in range(count_corners * 2):
                if j % 2 == 0:
                    x = int(corner_info[j])
                else:
                    y = int(corner_info[j])
                    corner_list.append(Point(x, y))
                    all_vertices.append(Point(x, y))

            Input.polygon_list.append(Polygon(corner_list, i + 1))

        start_information = lines[len(lines) - 1].split(" ")
        Input.start_node = Node(Point(int(start_information[0]), int(start_information[1])), -1)

        # Wenn Punkte doppelt vorkommen, haben Sie den Wert true, ansonsten false
        duplicates = {}
        for vertice in all_vertices:
            duplicates[vertice] = vertice in duplicates
        # Alle relevanten Eckpunkte, die nicht in mehr als einem Polygon vorkommen
        doulbe_corners = [i for i in duplicates if duplicates[i]]
        Input.double_corners.extend(doulbe_corners)


class Output:

    output_file = ""
    highest_time_after_departure = 0


    @staticmethod
    def display_result(result):
        final_output = "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' viewBox='0 0 1024 768'> " \
                       "<g transform='scale(1 -1)'> <g transform='translate(0 -768)'>" \
                       "<line id='y' x1='0' x2='0' y1='0' y2='768' fill='none' stroke='#212121' stroke-width='3'/>"
        final_output += Output.draw_polygons()

        if result is None:
            print ("Es gibt keinen Pfad zum Ziel")
            return

        shortest_path = result[0][0]
        visibility_graph = result[1]



        bus_departure_time = 27000      # 7:30 in Sekunden
        lisas_leaving = int(round(Output.highest_time_after_departure))     # Lisas Verlassen des Hauses
        start_time = str(datetime.timedelta(seconds=bus_departure_time+lisas_leaving))
        print ("Startzeit: " + start_time)


        all_nodes = visibility_graph.keys()

        output_lines = ""
        for node in all_nodes:
            if node.id == -1:
                # Kreis zeichnen fuer Start-Knoten
                output_lines += "<circle cx='" + str(node.point.x) +"' cy='" + \
                                str(node.point.y) +"' r='3' stroke='black' stroke-width='3' fill='red'></circle>"

            if node.id == -2:
                end_time = int(bus_departure_time + node.point.y*(3.0/25.0))
                print ("Zielzeit: " + str(datetime.timedelta(seconds=end_time)))      # Zielzeit
                print ("y-Koordinate des Auftreffens: " + str(node.point.y))        # y-Koordinaten des Auftreffens
                # Kreis zeichnen fuer End-Knoten
                output_lines += "<circle cx='" + str(node.point.x) +"' cy='" + \
                                str(node.point.y) +"' r='3' stroke='black' stroke-width='3' fill='blue'></circle>"
            for visible_node in visibility_graph[node]:
                output_lines += "\n <line "
                output_lines += "x1='" + str(node.point.x) + "' y1='" + str(node.point.y) + "'"
                output_lines += " x2='" + str(visible_node.point.x) + "' y2='" + str(visible_node.point.y) + "' "
                output_lines += "stroke='green' stroke-width='1'/>"

        length_lisas_path = int(round(result[0][1]))  # Laenge von Lisas Route in Minuten
        print ("Länge von Lisas Route: " + str(length_lisas_path) + " Meter")

        duration_lisas_path = round((length_lisas_path * (6.0 / 25.0)) / 60.0, 2)  # Dauer von Lisas Route in Minuten
        print ("Dauer von Lisas Route: " + str(duration_lisas_path) + " Minuten")

        lisas_route = ""
        for i, first_node in reversed(list(enumerate(shortest_path))):
            lisas_route += str(first_node) + " -->"
            if i > 0:
                second_node = shortest_path[i-1]
                output_lines += "\n <line "
                output_lines += "x1='" + str(first_node.point.x) + "' y1='" + str(first_node.point.y) + "'"
                output_lines += " x2='" + str(second_node.point.x) + "' y2='" + str(second_node.point.y) + "' "
                output_lines += "stroke='red' stroke-width='3'/>"

        print("Lisas Route: " + lisas_route + " Ziel erreicht")




        final_output += output_lines
        final_output += "</g> </g> </svg>"

        if Output.output_file in os.listdir("output"):
            os.remove("output/" + Output.output_file)

        with open("output/" + Output.output_file, "a") as f:
            f.write(final_output)
            f.close()

        print ("\nDie generierte SVG-Datei befindet sich im Output-Ordner")


    @staticmethod
    def draw_polygons():
        output = ""
        # Eckpunkte die zu zwei Polygonen gehoeren werden ausgefiltert
        all_vertices = []
        for polygon in Input.polygon_list:
            for vertice in polygon.points:
                all_vertices.append(vertice)

        for polygon in Input.polygon_list:
            output += "<polygon id='P" + str(polygon.id) + "' points='"
            for point in polygon.points:
                output += str(point.x) + " " + str(point.y) + " "
            output += "' fill='#6B6B6B' stroke='#212121' stroke-width='2'/> \n"


        return output


class Interface:

    @staticmethod
    def start():
        print ("Welche Umgebung soll eingelesen werden?")
        all_text_files = []
        for file in os.listdir("beispieldaten"):
            if file.endswith(".txt"):
                all_text_files.append(file)

        for i in range(len(all_text_files)):
            print ("[" + str(i+1) + "] " + all_text_files[i])

        choice = int(input("\nAuswahl: "))-1
        print ("\nRoute wird berechnet... \n")
        Input.polygon_liste_einlesen(all_text_files[choice])
        Output.output_file = all_text_files[choice].replace(".txt", "") + ".svg"
        Output.display_result(LisaRennt.get_optimal_result())



Interface.start()




