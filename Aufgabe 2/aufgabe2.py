import os
from math import sin, cos, pi, sqrt, acos, tan


# Repraesentiert einen Punkt in Form einer Koordinate oder einen Vektor mit dessen x- und y-Wert
class Point:
    def __init__(self, x, y):
        self.x = x      # x-Koordinate
        self.y = y      # y-Koordinate

    # Punkt als Zeichenkette
    def __repr__(self):
        return " (" + str(round(self.x, 2)) + " | " + str(round(self.y, 2)) + ") "

    # Vergleich-Funktion
    def __eq__(self, other):
        return (self.x, self.y) == (other.x, other.y)

    # Euklidische Distanz zweier Punkte berechnen (Anwenden des Satzes des Pythagoras)
    def euclidean_distance(self, other_point):
        cathetus_a = other_point.x - self.x
        cathetus_b = other_point.y - self.y
        hypotenuse = sqrt(cathetus_a ** 2 + cathetus_b ** 2)
        return hypotenuse

    # Multipliziert einen Ortsvektor mit einem bestimmten Faktor
    def scale(self, factor):
        return Point(self.x*factor, self.y*factor)


# Repraesentiert ein Dreieck mit einer Seite und dessen zwei anliegenden Winkel
class Triangle:
    def __init__(self, c, alpha, beta, id):
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.id = id

        self.point_a = None
        self.point_b = None
        self.point_c = None
        self.heuristic = 0

    def gamma(self):
        return pi - self.alpha - self.beta

    def construct_point_c(self):
        # 1. Laenge von BC mit Sinussatz berechnen
        length_bc = abs((self.c * sin(self.alpha)) / sin(self.gamma()))
        # 2. Zu rotierenden Vektor mit Laenge BC in Richtung c
        ratio = length_bc/self.c
        vector_c = Point(self.point_a.x - self.point_b.x, self.point_a.y - self.point_b.y)
        to_rotate = Point(ratio * vector_c.x, ratio * vector_c.y)
        angle = 2*pi - self.beta
        # 3. Ortsvektor OC mithilfe der Rotationsmatrix berechnen
        rotated_x = to_rotate.x * cos(angle) - to_rotate.y * sin(angle)
        rotated_y = to_rotate.x * sin(angle) + to_rotate.y * cos(angle)
        self.point_c = Point(rotated_x + self.point_b.x, rotated_y + self.point_b.y)

    def construct_point_a(self, target_point):
        # Erforderlich beim Beginn einer neuen Phase
        # Konstruktion des Punktes A, ausgehend vom Punkt B zum Zielpunkt
        ratio = self.c / self.point_b.euclidean_distance(target_point)
        prev_vector_bc = Point(target_point.x - self.point_b.x,
                               target_point.y - self.point_b.y)
        vector_ba = prev_vector_bc.scale(ratio)
        self.point_a = Point(self.point_b.x + vector_ba.x,
                                self.point_b.y + vector_ba.y)

    def __repr__(self):
        return str([self.point_a, self.point_b, self.point_c])

    # Vergleichsoperationen zum Sortieren nach den Heuristiken
    def __eq__(self, other):
        return self.heuristic == other.heuristic

    def __lt__(self, other):
        return self.heuristic > other.heuristic


# Repraesentiert eine Kante, bzw. eine Strecke mit zwei definierten Punkten A und B
class Edge:
    def __init__(self, point_a, point_b):
        self.point_a = point_a
        self.point_b = point_b

    # Prueft, ob zwei Kanten sich schneiden, wobei beide Kanten als Geraden fortgefuehrt werden
    def does_intersect_with(self, other_edge):
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
            ratio_edge_1 = self.get_ratio(Point(intersection_x, intersection_y))
            ratio_edge_2 = other_edge.get_ratio(Point(intersection_x, intersection_y))
            return 0 < ratio_edge_1 < 1 and 0 < ratio_edge_2 < 1

        # Steigung m1 kann berechnet werden, da dx1 ungleich 0 ist
        m1 = dy1 / dx1
        b1 = self.point_a.y - m1 * self.point_a.x

        # Wenn Kante b parallel zur y-Achse ist
        if dx2 == 0.0:
            intersection_x = other_edge.point_a.x
            intersection_y = m1 * intersection_x + b1
            ratio_edge_1 = self.get_ratio(Point(intersection_x, intersection_y))
            ratio_edge_2 = other_edge.get_ratio(Point(intersection_x, intersection_y))
            return 0 < ratio_edge_1 < 1 and 0 < ratio_edge_2 < 1

        # Steigung m2 kann berechnet werden, da dx2 ungleich 0 ist
        m2 = dy2 / dx2
        b2 = other_edge.point_a.y - m2 * other_edge.point_a.x
        intersection_x = (b2 - b1) / float(m1 - m2)
        intersection_y = m1 * intersection_x + b1
        ratio_edge_1 = self.get_ratio(Point(intersection_x, intersection_y))
        ratio_edge_2 = other_edge.get_ratio(Point(intersection_x, intersection_y))
        return 0 < ratio_edge_1 < 1 and 0 < ratio_edge_2 < 1

    # Erhalte den Faktor mit dem die Kante erweitert werden muss, um auf einen bestimmten Punkt zu treffen
    def get_ratio(self, point):
        # Zu pruefende Kante ist parallel zur Y-Achse --> Verhaeltnis ueber y-Werte
        if self.point_b.x == self.point_a.x:
            ratio = (point.y - self.point_a.y) / float(self.point_b.y - self.point_a.y)
        # Sonst: Verhaeltnis ueber x-Werte
        else:
            ratio = (point.x - self.point_a.x) / float(self.point_b.x - self.point_a.x)

        return ratio


# Klasse, die alle wichtigen Methoden zum Anordnen von Dreiecken enthaelt
class TrianglePlacement:
    remaining_angle = 0     # Verbleibende Winkel einer Phase
    new_phase = False       # True, wenn erstes Dreieck einer neuen

    @staticmethod
    def arrange_triangles(sorted_triangles):
        output = []

        # Anlegen des ersten Dreiecks
        current = sorted_triangles[0]
        sorted_triangles.pop(0)
        current.point_a = Point(0, 0)
        current.point_b = Point(current.c, 0)
        current.id = 1
        current.construct_point_c()
        TrianglePlacement.remaining_angle = pi - current.beta
        output.append(current)

        # Methode zum Setzen aller folgenden Dreiecke
        return TrianglePlacement.place_triangle(current, sorted_triangles, output)

    # Setzen eines Dreiecks an die Kante des vorherigen
    @staticmethod
    def place_triangle(previous, to_place, output):
        # Wenn ein Dreieck 90 Grad der Phase ueberschreitet...
        if TrianglePlacement.remaining_angle - to_place[0].beta <= pi/2:
            # ... Liste von hinten abarbeiten
            current = to_place[len(to_place)-1]
            to_place.pop(len(to_place)-1)

        else:
            # Ansonsten von vorne abarbeiten
            current = to_place[0]
            to_place.pop(0)

        current.id = previous.id + 1    # Setzen der ID

        if current.beta > TrianglePlacement.remaining_angle:
            # Eine neue Phase muss gestartet werden
            # Das Dreieck wird mit c so nah wie moeglich ans vorherige geschoben
            to_place.append(current)
            current = to_place[len(to_place)-1]
            to_place.pop(len(to_place)-1)

            # Tausche Beta mit alpha
            save_alpha = current.alpha
            current.alpha = current.beta
            current.beta = save_alpha

            # Zu setzende Kante als Gerade: y(x) = m * x + b
            # y(x) = tan^-1(beta) * x + (c_y - m * c_x)
            m = tan(current.alpha)
            b = previous.point_c.y - m * previous.point_c.x

            # x-Koordinate von Punkt A ist Nullstelle von y(x)
            # 0 = m*x +b <=> x_n = -b/m
            point_a_x = -b / m

            current.point_a = Point(point_a_x, 0)
            current.point_b = Point(current.point_a.x + current.c, 0)
            current.construct_point_c()
            output.append(current)
            TrianglePlacement.remaining_angle = pi - current.beta
            TrianglePlacement.new_phase = True  # Naechstes Dreieck kann andere Dreiecke schneiden
            Output.last_point = current.point_a

        else:
            current.point_b = previous.point_b
            if TrianglePlacement.new_phase:
                # Beim ersten Dreieck einer neuen Phase (sofern es nicht die erste ist)
                # muss das Dreieck eventuell justiert werden, damit es nicht
                # zu Ueberschneidungen kommt
                TrianglePlacement.new_phase = False
                current.construct_point_a(previous.point_c)
                additional_angle = 0
                # Iteriere durch jedes gesetzte Dreieck, bis auf das letzte
                for triangle in output[:len(output)-1]:
                    # Aktuelle Kante c
                    edge_c_new = Edge(current.point_b, current.point_a)
                    prev_edge_a = Edge(triangle.point_b, triangle.point_c)  # Vorherige Kante a
                    prev_edge_b = Edge(triangle.point_c, triangle.point_a)  # Vorherige Kante b

                    if prev_edge_a.does_intersect_with(edge_c_new):
                        # Neue Kante c schneidet vorherige Kante a --> Justierung erforderlich
                        current.construct_point_a(triangle.point_c)

                        # Winkel zwischen BC und BAnew berechnen
                        vector_bc = Point(previous.point_c.x - previous.point_b.x,
                                          previous.point_c.y - previous.point_b.y)
                        vector_ba_new = Point(current.point_a.x - previous.point_b.x,
                                          current.point_a.y - previous.point_b.y)
                        scalar = (vector_ba_new.x * vector_bc.x + vector_ba_new.y * vector_bc.y)
                        length_product = sqrt(vector_ba_new.y**2 + vector_ba_new.x**2) * \
                                         sqrt(vector_bc.y**2 + vector_bc.x**2)
                        additional_angle = acos(scalar/length_product)

                    edge_c_new = Edge(current.point_b, current.point_a)

                    if prev_edge_b.does_intersect_with(edge_c_new):
                        # Neue Kante c schneidet vorherige Kante b --> Justierung erforderlich
                        current.construct_point_a(triangle.point_a)

                        # Winkel zwischen BC und BAnew berechnen
                        vector_bc = Point(previous.point_c.x - previous.point_b.x,
                                          previous.point_c.y - previous.point_b.y)
                        vector_ba_new = Point(current.point_a.x - previous.point_b.x,
                                              current.point_a.y - previous.point_b.y)
                        scalar = (vector_ba_new.x * vector_bc.x + vector_ba_new.y * vector_bc.y)
                        length_product = sqrt(vector_ba_new.y ** 2 + vector_ba_new.x ** 2) * \
                                         sqrt(vector_bc.y ** 2 + vector_bc.x ** 2)

                        additional_angle = acos(scalar / length_product)

                # Verbleibenden Winkel aktualisieren
                TrianglePlacement.remaining_angle -= current.beta + additional_angle
            else:
                TrianglePlacement.remaining_angle -= current.beta
                # Berechne Vektor von B zu A (Kante c) durch Skalieren der Seite BC des vorherigen Dreiecks
                current.construct_point_a(previous.point_c)

            # Konstruktion des neuen Punktes C
            current.construct_point_c()
            output.append(current)
            Output.last_point = current.point_b

        if len(to_place) == 0:
            return output

        return TrianglePlacement.place_triangle(current, to_place, output)


class Input:
    triangle_list = []

    @staticmethod
    def read_triangles(file_name):
        lines = [line.rstrip('\n') for line in open('beispieldaten/' + file_name)]
        count_triangles = int(lines[0])
        lines.pop(0)
        all_vertices = []
        total_length = 0        # Gesamtlaenge aller laengsten Seiten eines Dreiecks
        angle_sum = 0           # Summe aller kleinsten Winkel
        for i in range(count_triangles):
            corner_info = lines[i].split(" ")
            corner_info.pop(0)
            corner_list = []
            x = 0
            for j in range(6):
                if j % 2 == 0:
                    x = int(corner_info[j])
                else:
                    y = int(corner_info[j])
                    corner_list.append(Point(x, y))
                    all_vertices.append(Point(x, y))

            # Seite c des Dreiecks ermitteln (laengste Seite) + Punkte hinzufuegen
            side_a = corner_list[0].euclidean_distance(corner_list[1])
            side_b = corner_list[0].euclidean_distance(corner_list[2])
            side_c = corner_list[1].euclidean_distance(corner_list[2])
            current_sides = [side_a, side_b, side_c]

            biggest_side = max(current_sides)

            if side_a == biggest_side:
                current_sides[0] = current_sides[2]
                current_sides[2] = biggest_side
            if side_b == biggest_side:
                current_sides[1] = current_sides[2]
                current_sides[2] = biggest_side

            total_length += biggest_side     # Zur Berechnung der Heuristik

            # Berechnen von Alpha und Beta (der kleinere Winkel der beiden)
            # Anwenden des Kosinussatz
            side_b = current_sides[1]
            side_a = current_sides[0]
            side_c = current_sides[2]
            alpha = acos((side_b**2 + side_c**2 - side_a**2) / (2*side_b*side_c))
            beta = acos((side_c**2 + side_a**2 - side_b**2) / (2*side_c*side_a))

            # Beta soll der kleinere Winkel sein
            if alpha < beta:
                copy_alpha = alpha
                alpha = beta
                beta = copy_alpha

            angle_sum += beta

            current_triangle = Triangle(biggest_side, alpha, beta, i)
            Input.triangle_list.append(current_triangle)

        # Zuweisen der Heuristiken
        for triangle in Input.triangle_list:
            quota_length = triangle.c / total_length
            quota_angle = triangle.beta / angle_sum
            triangle.heuristic = quota_length - quota_angle

        # Sortieren der Dreiecke nach Heuristik
        Input.triangle_list = sorted(Input.triangle_list)


class Output:
    output_file = "test.svg"
    last_point = None

    @staticmethod
    def display_result(result):
        final_output = "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' viewBox='0 0 1024 768'> \n" \
                       "<g transform='scale(1 -1)'> <g transform='translate(0 -600)' fill='#ffcc99'>\n" \
                       "<line id='y' x1='0' x2='820' y1='0' y2='0' fill='none' stroke='#000000' stroke-width='3'/>"
        if Output.last_point is not None:
            gesamtabstand = Output.last_point.x - result[0].point_b.x
        else:
            gesamtabstand = 0
        print ("Gesamtsbtand " + str(round(gesamtabstand, 2)))
        print("Anordnung der Dreiecke:")
        for triangle in result:
            final_output += Output.draw_triangle(triangle)
            print("D" + str(triangle.id) + ": " + "A" + str(triangle.point_a) + ", B" + str(triangle.point_b)
                  + ", C" + str(triangle.point_c))
        final_output += "</g> </g> </svg>"

        if Output.output_file in os.listdir("output"):
            os.remove("output/" + Output.output_file)

        with open("output/" + Output.output_file, "a") as f:
            f.write(final_output)
            f.close()

        print ("\nDie generierte SVG-Datei befindet sich im Output-Ordner \n ")

    @staticmethod
    def draw_triangle(triangle):
        output = ""
        output += "<polygon id='D" + str(triangle.id) + "' points='"
        for vertice in [triangle.point_a, triangle.point_b, triangle.point_c]:
            output += str(round(vertice.x, 2)) + " " + str(round(vertice.y, 2)) + " "
        output += "' stroke='#212121' stroke-width='2'/> \n"

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
        Input.read_triangles(all_text_files[choice])
        Output.output_file = all_text_files[choice].replace(".txt", "") + ".svg"


Interface.start()
sorted_triangles = Input.triangle_list
result = TrianglePlacement.arrange_triangles(sorted_triangles)
Output.display_result(result)



