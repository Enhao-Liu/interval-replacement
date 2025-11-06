import sympy as sp
from itertools import combinations

# --- Helper functions ---
def powerset(elements):
    """Generate all possible subsets of a set of elements."""
    all_subsets = []
    for subset_size in range(len(elements) + 1):
        all_subsets.extend(combinations(elements, subset_size))
    return [list(subset) for subset in all_subsets]

def flatten(xss):
    """Flatten a list of lists."""
    return [x for xs in xss for x in xs]

# --- Interval Class ---
class Interval:
    """Define intervals by their sources and sinks."""
    def __init__(self, src, snk):
        self.src = [tuple(s) for s in src]
        self.snk = [tuple(s) for s in snk]

# --- Representation Class ---
class Representation:
    """Quiver representation with symbolic matrices using SymPy."""
    def __init__(self, dimensions):
        self.dimensions = dimensions
        self.nodes = self.generate_nodes()
        self.edges = self.generate_edges()
        self.vecs = {}  # dim of vector space at node
        self.mats = {}  # linear maps (symbolic matrices)

    # --- Grid and edges ---
    def generate_nodes(self):
        nodes = []
        for coordinates in self.generate_coordinates():
            nodes.append(tuple(coordinates))
        return nodes

    def generate_coordinates(self, current_dim=0, current_coords=None):
        if current_coords is None:
            current_coords = []
        if current_dim == len(self.dimensions):
            return [tuple(current_coords)]
        coordinates = []
        for i in range(self.dimensions[current_dim]):
            new_coords = current_coords + [i]
            coordinates.extend(self.generate_coordinates(current_dim + 1, new_coords))
        return coordinates

    def generate_edges(self):
        edges = []
        for node1 in self.nodes:
            for node2 in self.nodes:
                if all(node1[i] <= node2[i] for i in range(len(self.dimensions))):
                    edges.append((node1, node2))
        return edges

    # --- Vector spaces and matrices ---
    def create_vecs(self, node, m):
        if node in self.nodes:
            self.vecs[node] = m
        else:
            print("Error: node not in grid.")

    def create_matrix(self, node1, node2, matrix):
        if (node1, node2) in self.edges:
            u = self.vecs.get(node1, 0)
            v = self.vecs.get(node2, 0)
            if matrix is None:
                # store a zero matrix of the right size
                matrix = sp.zeros(v, u)
            else:
                if matrix.shape != (v, u):
                    print(f"Error: matrix dimensions mismatch for edge {node1}->{node2}. Expected ({v},{u}), got {matrix.shape}.")
                    return
            self.mats[(node1, node2)] = matrix
        else:
            print(f"Error: nodes {node1} and {node2} are not connected.")

    # --- Lattice operations ---
    def join(self, node1, node2):
        return tuple(max(a, b) for a, b in zip(node1, node2))

    def meet(self, node1, node2):
        return tuple(min(a, b) for a, b in zip(node1, node2))

    def is_smaller(self, node1, node2):
        return all(node1[i] <= node2[i] for i in range(len(self.dimensions)))

    # --- Convex hull of an interval ---
    def convex_square(self, node1, node2):
        if not self.is_smaller(node1, node2):
            return []
        points = []
        for coordinates in self.nodes:
            if all(node1[i] <= coordinates[i] <= node2[i] for i in range(len(self.dimensions))):
                points.append(tuple(coordinates))
        return points

    def int_hull(self, interval):
        points = []
        for src in interval.src:
            for snk in interval.snk:
                points += self.convex_square(src, snk)
        points = list(set(points))
        for pt in points[:]:
            if self.to_be_removed(pt, interval):
                points.remove(pt)
        return points

    def to_be_removed(self, point, interval):
        for src in interval.src:
            if self.is_smaller(point, src) and point != src:
                return True
        for snk in interval.snk:
            if self.is_smaller(snk, point) and point != snk:
                return True
        return False

    # --- Sources and sinks ---
    def get_src_snk(self, points):
        src, snk = [], []
        for pt in points:
            if self.add_source(pt, points):
                src.append(pt)
            if self.add_sink(pt, points):
                snk.append(pt)
        return src, snk

    def add_source(self, pt, points):
        for other in points:
            if self.is_smaller(other, pt) and other != pt:
                return False
        return True

    def add_sink(self, pt, points):
        for other in points:
            if self.is_smaller(pt, other) and other != pt:
                return False
        return True

    # --- Evaluation of a path ---
    def evaluation(self, node1, node2):
        # handle empty vector spaces
        dim_out = self.vecs.get(node2, 0)
        dim_in = self.vecs.get(node1, 0)
        if dim_out == 0 or dim_in == 0:
            return sp.zeros(dim_out, dim_in)

        # start with identity matrix of appropriate size
        mat = sp.eye(dim_out)
        current = list(node2)
        previous = list(node2)

        for i in range(len(self.dimensions)):
            k = current[i] - node1[i]
            while k > 0:
                previous[i] -= 1
                edge = (tuple(previous), tuple(current))
                if self.mats.get(edge) is not None:
                    mat = mat * self.mats[edge]  # symbolic multiplication
                else:
                    return sp.zeros(dim_out, dim_in)
                current[i] -= 1
                k -= 1

        return mat


    # --- Construct block matrices ---
    def construct_matrix_MN(self, column_labels, block_signature, dual = False):
        '''
        Construct a matrix according to the block_signature: a list with elements [row, [col1, col2]].
        For each we add a block row with M_{col1, row} and -M_{col2, row}. The parameter column_labels
        determines the order of the column blocks.
        If dual is True, we add M_{row, col1}^T and -M_{row, col2}^T, and return transpose.
        '''
        col_num = sum(self.vecs[node] for node in column_labels)
        row_num = sum(self.vecs[node] for node, _ in block_signature)

        idx_col = {}  # column indices where each block should start and end
        start = 0
        for node in column_labels:
            idx_col[node] = (start, start + self.vecs[node])
            start = idx_col[node][1]

        matrix = sp.zeros(row_num, col_num)
        r_current, r_next = 0, 0  # current row

        for row, columns in block_signature:
            r_next = r_current + self.vecs[row]
            for i, col in enumerate(columns):
                block = self.evaluation(col, row) if not dual else self.evaluation(row, col).T
                matrix[r_current : r_next, idx_col[col][0] : idx_col[col][1]] = (-1)**i * block
            r_current = r_next

        return matrix if not dual else matrix.T


    def construct_matrix_M_tot(self, interval):
        '''Given an interval with n sources, return the matrix_M.'''
        column_labels = interval.src
        block_signature = []
        for i, src1 in enumerate(interval.src):
            for src2 in interval.src[i + 1:]:
                block_signature.append((self.join(src1, src2), (src1, src2)))
        return self.construct_matrix_MN(column_labels, block_signature)


    def construct_matrix_N_tot(self, interval):
        '''Given an interval with n sinks, return the matrix_N.'''
        column_labels = interval.snk
        block_signature = []
        for i, snk1 in enumerate(interval.snk):
            for snk2 in interval.snk[i + 1:]:
                block_signature.append((self.meet(snk1, snk2), (snk1, snk2)))
        return self.construct_matrix_MN(column_labels, block_signature, dual=True)


    def construct_matrix_M_ss(self, interval):
        block_signature = [(bound, (src1, src2))
                           for i, src1 in enumerate(interval.src)
                           for src2 in interval.src[i+1:]
                           for bound in self.find_upper_bounds(interval, src1, src2)]
        return self.construct_matrix_MN(interval.src, block_signature)

    def construct_matrix_N_ss(self, interval):
        block_signature = [(bound, (snk1, snk2))
                           for i, snk1 in enumerate(interval.snk)
                           for snk2 in interval.snk[i+1:]
                           for bound in self.find_lower_bounds(interval, snk1, snk2)]
        return self.construct_matrix_MN(interval.snk, block_signature, dual=True)

    def find_source_sink_indices_with_path(self, interval):
        for i in range(len(interval.src)):
            for j in range(len(interval.snk)):
                if self.is_smaller(interval.src[i], interval.snk[j]):
                    return i, j
        raise ValueError("No source -> sink path found.")
    
    def find_upper_bounds(self, interval, node1, node2):
        """
        Return a list of all sinks of the given interval larger than both given nodes.
        """
        return [snk for snk in interval.snk if self.is_smaller(node1, snk) and self.is_smaller(node2, snk)]

    def find_lower_bounds(self, interval, node1, node2):
        """
        Return a list of all sources of the given interval smaller than both given nodes.
        """
        return [src for src in interval.src if self.is_smaller(src, node1) and self.is_smaller(src, node2)]

    def list_int(self, conv=True):
        '''
        Return the list of all intervals possible.

        If conv == True, return intervals in the form of convex hulls.
        If conv == False, return intervals in the form Interval(src, snk).
        '''
        nodes = self.nodes
        candidates = powerset(nodes)[1:]  # remove the empty interval
        list_int = []
        for c in candidates:
            if self.is_convex(c) and self.is_connected(c):
                if conv:
                    list_int.append(c)
                else:
                    tmp = self.get_src_snk(c)
                    list_int.append(Interval(tmp[0], tmp[1]))
        return list_int

    def int_rank(self, interval, compression='tot'):
        i, j = self.find_source_sink_indices_with_path(interval)
        new_src, new_snk = interval.src.copy(), interval.snk.copy()
        new_src[0], new_src[i] = new_src[i], new_src[0]
        new_snk[0], new_snk[j] = new_snk[j], new_snk[0]
        interval = Interval(new_src, new_snk)

        if compression == 'tot':
            M = self.construct_matrix_M_tot(interval)
            N = self.construct_matrix_N_tot(interval)
            mat = self.evaluation(interval.src[0], interval.snk[0])
        elif compression == 'ss':
            M = self.construct_matrix_M_ss(interval)
            N = self.construct_matrix_N_ss(interval)
            mat = self.evaluation(interval.src[0], interval.snk[0])
        else:
            raise ValueError("Invalid compression")

        M_rows, M_cols = M.rows if M.rows > 0 else 1, M.cols if M.cols > 0 else 1
        N_rows, N_cols = N.rows if N.rows > 0 else 1, N.cols if N.cols > 0 else 1

        # Safe C matrix (match mat shape)
        C = sp.zeros(max(mat.rows, N_rows), max(mat.cols, M_cols))
        C[:mat.rows, :mat.cols] = mat

        # Safe zero blocks
        M_safe = M if M.rows > 0 and M.cols > 0 else sp.zeros(M_rows, M_cols)
        N_safe = N if N.rows > 0 and N.cols > 0 else sp.zeros(N_rows, N_cols)
        zeros_top_right = sp.zeros(M_rows, N_cols)

        # Construct block matrix safely
        block = sp.BlockMatrix([[M_safe, zeros_top_right], [C, N_safe]]).as_explicit()

        def symbolic_rank(matrix):
            """
            Compute the symbolic rank of a (possibly rectangular) SymPy matrix.
            Returns a sympy expression (Piecewise) depending on symbolic parameters.
            """
            rows, cols = matrix.shape
            if rows == 0 or cols == 0:
                return 0

            # The rank is at most min(rows, cols)
            max_rank = min(rows, cols)

            # Check all square minors from largest to smallest
            for k in range(max_rank, 0, -1):
                # Generate all sets of k rows and k columns
                from itertools import combinations
                row_combinations = list(combinations(range(rows), k))
                col_combinations = list(combinations(range(cols), k))

                for r_idx in row_combinations:
                    for c_idx in col_combinations:
                        minor = matrix.extract(list(r_idx), list(c_idx))
                        det = minor.det().simplify()
                        if not det.is_zero:
                            return sp.Piecewise((k, sp.Ne(det, 0)), (0, True))
            return 0

        return symbolic_rank(block) - symbolic_rank(M_safe) - symbolic_rank(N_safe)

    # --- Interval replacement ---
    def int_replacement(self, interval, compression='tot'):
        repl = self.int_rank(interval, compression)
        cov_ps = powerset(self.cover(interval))
        for c in cov_ps[1:]:
            if len(c) == 1:
                tmp = self.get_src_snk(c[0])
                i = Interval(tmp[0], tmp[1])
                repl -= self.int_rank(i, compression)
            if len(c) > 1:
                eps = len(c)
                c_flat = [list(set(flatten(c)))]
                tmp = self.get_src_snk(c_flat[0])
                i2 = Interval(tmp[0], tmp[1])
                repl += (-1)**eps * self.int_rank(i2, compression)
        return(repl)

    def cover(self, interval, conv=True):
        """Compute the cover of an interval without duplicates."""
        cover_list = []
        hull = self.int_hull(interval)
        seen = set()  # Track seen intervals to avoid duplicates
        
        for point in hull:
            if self.is_corner(point, interval):
                # Try extending in each dimension
                for i in range(len(self.dimensions)):
                    # Try increasing coordinate
                    if point[i] < self.dimensions[i] - 1:
                        new_point = list(point)
                        new_point[i] += 1
                        new_point = tuple(new_point)
                        
                        if new_point not in hull:
                            # Create new hull by adding the new point
                            new_hull = hull + [new_point]
                            if self.is_convex(new_hull):
                                if conv:
                                    # For point lists: use sorted tuple as signature
                                    signature = tuple(sorted(new_hull))
                                    if signature not in seen:
                                        cover_list.append(new_hull)
                                        seen.add(signature)
                                else:
                                    # For Interval objects: get sources/sinks
                                    tmp = self.get_src_snk(new_hull)
                                    new_interval = Interval(tmp[0], tmp[1])
                                    sig = (tuple(sorted(new_interval.src)), tuple(sorted(new_interval.snk)))
                                    if sig not in seen:
                                        cover_list.append(new_interval)
                                        seen.add(sig)
                    
                    # Try decreasing coordinate  
                    if point[i] > 0:
                        new_point = list(point)
                        new_point[i] -= 1
                        new_point = tuple(new_point)
                        
                        if new_point not in hull:
                            # Create new hull by adding the new point
                            new_hull = hull + [new_point]
                            if self.is_convex(new_hull):
                                if conv:
                                    # For point lists: use sorted tuple as signature
                                    signature = tuple(sorted(new_hull))
                                    if signature not in seen:
                                        cover_list.append(new_hull)
                                        seen.add(signature)
                                else:
                                    # For Interval objects: get sources/sinks
                                    tmp = self.get_src_snk(new_hull)
                                    new_interval = Interval(tmp[0], tmp[1])
                                    sig = (tuple(sorted(new_interval.src)), tuple(sorted(new_interval.snk)))
                                    if sig not in seen:
                                        cover_list.append(new_interval)
                                        seen.add(sig)
        
        return cover_list

    def is_convex(self, points):
        tmp = self.get_src_snk(points)
        interval = Interval(tmp[0], tmp[1])
        return set(points) == set(self.int_hull(interval))

    def is_corner(self, point, interval):
        for i in range(len(self.dimensions)):
            aligned = False
            for src in interval.src:
                if src[i] == point[i]:
                    aligned = True
            for snk in interval.snk:
                if snk[i] == point[i]:
                    aligned = True
            if not aligned:
                return False
        return True

    def is_connected(self, points):
        visited = set()
        def dfs(node):
            visited.add(node)
            for neighbor in points:
                if neighbor not in visited and (self.is_smaller(node, neighbor) or self.is_smaller(neighbor, node)):
                    dfs(neighbor)
        dfs(points[0])
        return all(node in visited for node in points)
    
    # --- Check validity of symbolic calculation ---
    def find_complete_singularities(self, parameter):
        """
        Find all singular parameter values where interval replacement computations may change.
        This identifies where symbolic computation cannot be simply substituted.
        """
        print("Identifying singular parameter values for interval replacement")
        
        all_intervals = self.list_int(conv=False)
        print(f"Total intervals in quiver: {len(all_intervals)}")
        
        all_replacement_intervals = set()
        
        for interval in all_intervals:
            all_replacement_intervals.add((tuple(sorted(interval.src)), tuple(sorted(interval.snk))))
            
            cover_intervals = self.cover(interval, conv=False)
            for cov in cover_intervals:
                all_replacement_intervals.add((tuple(sorted(cov.src)), tuple(sorted(cov.snk))))
            
            cov_ps = powerset(cover_intervals)
            for c in cov_ps:
                if len(c) > 1:
                    all_points = []
                    for interval_obj in c:
                        hull_points = self.int_hull(interval_obj)
                        all_points.extend(hull_points)
                    unique_points = list(set(all_points))
                    tmp = self.get_src_snk(unique_points)
                    combined_int = Interval(tmp[0], tmp[1])
                    all_replacement_intervals.add((tuple(sorted(combined_int.src)), tuple(sorted(combined_int.snk))))
        
        print(f"Total unique intervals in replacement formulas: {len(all_replacement_intervals)}")
        
        complete_singularities = set()
        
        def construct_block_matrix(rep, src, snk):
            interval = Interval(list(src), list(snk))
            i, j = rep.find_source_sink_indices_with_path(interval)
            new_src, new_snk = interval.src.copy(), interval.snk.copy()
            new_src[0], new_src[i] = new_src[i], new_src[0]
            new_snk[0], new_snk[j] = new_snk[j], new_snk[0]
            interval = Interval(new_src, new_snk)

            M = rep.construct_matrix_M_ss(interval)
            N = rep.construct_matrix_N_ss(interval)
            mat = rep.evaluation(interval.src[0], interval.snk[0])

            M_rows, M_cols = M.rows if M.rows > 0 else 1, M.cols if M.cols > 0 else 1
            N_rows, N_cols = N.rows if N.rows > 0 else 1, N.cols if N.cols > 0 else 1

            C = sp.zeros(max(mat.rows, N_rows), max(mat.cols, M_cols))
            C[:mat.rows, :mat.cols] = mat

            M_safe = M if M.rows > 0 and M.cols > 0 else sp.zeros(M_rows, M_cols)
            N_safe = N if N.rows > 0 and N.cols > 0 else sp.zeros(N_rows, N_cols)
            zeros_top_right = sp.zeros(M_rows, N_cols)

            block = sp.BlockMatrix([[M_safe, zeros_top_right], [C, N_safe]]).as_explicit()
            return block
        
        def find_singularities(matrix, param):
            rows, cols = matrix.shape
            max_rank = min(rows, cols)
            singularities = set()
            
            for k in range(max_rank, 0, -1):
                for row_indices in combinations(range(rows), k):
                    for col_indices in combinations(range(cols), k):
                        minor = matrix.extract(list(row_indices), list(col_indices))
                        det = minor.det().simplify()
                        
                        if det != 0:
                            solutions = sp.solve(det, param)
                            for sol in solutions:
                                try:
                                    singularities.add(float(sol.evalf()))
                                except:
                                    singularities.add(sol)
                            return singularities
            return singularities
        
        for i, (src, snk) in enumerate(all_replacement_intervals):
            if i % 20 == 19:
                print(f"Progress: {i+1}/{len(all_replacement_intervals)}")
                
            try:
                block = construct_block_matrix(self, src, snk)
                singularities = find_singularities(block, parameter)
                complete_singularities.update(singularities)
            except:
                continue
        
        print(f"Singular parameter values: {complete_singularities}")
        
        return complete_singularities

    # --- Utility ---
    def elements(self):
        print("Nodes:", self.nodes)
        print("Edges:", self.edges)
        print("Vector spaces:", self.vecs)
        print("Matrices:", self.mats)
