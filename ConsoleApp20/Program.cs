using System;
using System.Collections;
using System.Collections.Generic;

public class IntGraphNode
{
    public int value;
    public List<IntGraphNode> neighbors;

    public IntGraphNode(int value = 0, List<IntGraphNode> neighbors = null)
    {
        this.value = value;
        this.neighbors = neighbors ?? new List<IntGraphNode>();
    }
}
//העתקת הגרף המקורי

public class Solution
{
    public Dictionary<int, List<int>> CopyGraph(IntGraphNode node)
    {
        // מילון לשמירת הצמתים שהועתקו
        var copied = new Dictionary<IntGraphNode, IntGraphNode>();

        // פונקציה חזורית שמבצעת חיפוש DFS
        IntGraphNode DFS(IntGraphNode currentNode)
        {
            if (currentNode == null) return null;

            if (copied.ContainsKey(currentNode))
                return copied[currentNode];

            // יצירת עותק חדש לצומת
            var copy = new IntGraphNode(currentNode.value);
            copied[currentNode] = copy;

            // מעתיקים את השכנים
            foreach (var neighbor in currentNode.neighbors)
            {
                copy.neighbors.Add(DFS(neighbor));
            }

            return copy;
        }

        // התחלת החיפוש מהצומת הנתון
        DFS(node);

        // יצירת המילון המייצג את הגרף
        var result = new Dictionary<int, List<int>>();
        foreach (var kvp in copied)
        {
            var neighborsList = new List<int>();
            foreach (var neighbor in kvp.Value.neighbors)
            {
                neighborsList.Add(neighbor.value);
            }
            result[kvp.Value.value] = neighborsList;
        }

        return result;
    }
}
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//לבצע חיפוש DFS (Depth First Search) כדי לבדוק אם הגרף מחובר. אם הצלחנו להגיע לכל הצמתים ללא פיצול, אז הגרף מחובר
public class Solution2
{
    public bool ValidTree(int n, int[][] edges)
    {
        // אם יש פחות מ-1 צומת או יותר מדי קשתות
        if (n - 1 != edges.Length)
        {
            return false;
        }

        // יצירת רשימת שכנים לכל צומת
        var graph = new Dictionary<int, List<int>>();

        // יצירת הגרף (המילון של הצמתים ושכנותיהם)
        foreach (var edge in edges)
        {
            if (!graph.ContainsKey(edge[0]))
                graph[edge[0]] = new List<int>();
            if (!graph.ContainsKey(edge[1]))
                graph[edge[1]] = new List<int>();

            graph[edge[0]].Add(edge[1]);
            graph[edge[1]].Add(edge[0]);
        }

        // מערך למעקב אחרי הצמתים שביקרנו בהם
        var visited = new bool[n];

        // פונקציית DFS
        bool DFS(int node, int parent)
        {
            visited[node] = true;

            // חיפוש בDFS על כל השכנים
            foreach (var neighbor in graph[node])
            {
                // אם השכן הוא ההורה של הצומת הנוכחי, דלג עליו
                if (neighbor == parent)
                    continue;

                // אם כבר ביקרנו בשכן זה, זאת אומרת שיש cycle
                if (visited[neighbor])
                    return false;

                // המשך בDFS
                if (!DFS(neighbor, node))
                    return false;
            }
            return true;
        }

        // בדוק אם אפשר להתחיל חיפוש מ-צומת 0
        if (!DFS(0, -1))
            return false;

        // אם לא כל הצמתים בוקרו, הגרף לא מחובר
        foreach (var v in visited)
        {
            if (!v)
                return false;
        }

        return true;
    }
}


/// ////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
        //חיפוש במטריציה
        // מהנקודה למעלה למטה ימינה שמולה
public class MatrixDFS
{
    // מטריצה לדוגמה
    private int[,] matrix;
    private bool[,] visited;

    public MatrixDFS(int[,] grid)
    {
        matrix = grid;
        visited = new bool[matrix.GetLength(0), matrix.GetLength(1)];
    }

    public void DFS(int row, int col)
    {
        // אם מחוץ לגבולות המטריצה
        if (row < 0 || col < 0 || row >= matrix.GetLength(0) || col >= matrix.GetLength(1))
            return;

        // אם כבר ביקרנו בתא או שהתא אינו מתאים
        if (visited[row, col] || matrix[row, col] == 0)
            return;

        // סימון התא כנבקר
        visited[row, col] = true;

        Console.WriteLine($"Visiting cell ({row}, {col}) with value {matrix[row, col]}");

        // חיפוש DFS בארבעת הכיוונים (למעלה, למטה, שמאלה, ימינה)
        DFS(row - 1, col); // למעלה
        DFS(row + 1, col); // למטה
        DFS(row, col - 1); // שמאלה
        DFS(row, col + 1); // ימינה
    }
}












public class Program
{
    public static void Main()
    {
        // יצירת צמתים לדוגמה
        IntGraphNode n1 = new IntGraphNode(1);
        IntGraphNode n2 = new IntGraphNode(2);
        IntGraphNode n3 = new IntGraphNode(3);
        IntGraphNode n4 = new IntGraphNode(4);

        n1.neighbors.Add(n2);
        n1.neighbors.Add(n4);
        n2.neighbors.Add(n1);
        n2.neighbors.Add(n3);
        n3.neighbors.Add(n2);
        n3.neighbors.Add(n4);
        n4.neighbors.Add(n1);
        n4.neighbors.Add(n3);

        // העתקת הגרף
        Solution solution = new Solution();
        var copiedGraph = solution.CopyGraph(n1);

        // הצגת התוצאה
        foreach (var kvp in copiedGraph)
        {
            Console.Write(kvp.Key + ": ");
            foreach (var neighbor in kvp.Value)
            {
                Console.Write(neighbor + " ");
            }
            Console.WriteLine();
        }
        ///////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////

        //לבצע חיפוש DFS (Depth First Search) כדי לבדוק אם הגרף מחובר. אם הצלחנו להגיע לכל הצמתים ללא פיצול, אז הגרף מחובר


        Solution2 solution2 = new Solution2();

        // דוגמת קלט 1
        int[][] edges1 = new int[][]
        {
            new int[] {0, 1},
            new int[] {2, 3}
        };
        Console.WriteLine(solution2.ValidTree(4, edges1));  // Output: false

        // דוגמת קלט 2 (עץ תקני)
        int[][] edges2 = new int[][]
        {
            new int[] {0, 1},
            new int[] {0, 2},
            new int[] {2, 3}
        };
        Console.WriteLine(solution2.ValidTree(4, edges2));  // Output: true



        ///////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        //חיפוש במטריציה
        // מהנקודה למעלה למטה ימינה שמולה


        // מטריצה לדוגמה
        int[,] grid = new int[,]
        {
            { 1, 0, 1 },
            { 1, 0, 0 },
            { 0, 0, 1 }
        };

        MatrixDFS dfs = new MatrixDFS(grid);

        Console.WriteLine("Starting DFS from (0, 0):");
        dfs.DFS(0, 0);

        Console.WriteLine("\nStarting DFS from (2, 2):");
        dfs.DFS(2, 2);
        ///////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////
        ///שינוי הצבע 0 לצבע 2 במטריצה במאצעות חיפוש DFS




        int[][] image = new int[][]
               {
            new int[] { 1, 0, 1 },
            new int[] { 1, 0, 0 },
            new int[] { 0, 0, 1 }
               };

        int sr = 1, sc = 1, newColor = 2;

        Console.WriteLine("Original Image:");
        PrintImage(image);

        int[][] result = FloodFill(image, sr, sc, newColor);

        Console.WriteLine("\nFlood-Filled Image:");
        PrintImage(result);



    }


    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///שינוי הצבע 0 לצבע 2 במטריצה במאצעות חיפוש DFS


    static int[][] FloodFill(int[][] image, int sr, int sc, int newColor)
    {
        int startColor = image[sr][sc];

        // אם הצבע ההתחלתי כבר שווה לצבע החדש, אין צורך לשנות
        if (startColor == newColor)
            return image;

        // קריאה לפונקציה רקורסיבית לביצוע Flood Fill
        DFS(image, sr, sc, startColor, newColor);

        return image;
    }

    static void DFS(int[][] image, int row, int col, int startColor, int newColor)
    {
        // בדיקת גבולות המטריצה
        if (row < 0 || row >= image.Length || col < 0 || col >= image[0].Length)
            return;

        // בדיקה אם התא הנוכחי הוא בצבע ההתחלתי
        if (image[row][col] != startColor)
            return;

        // צביעה בצבע החדש
        image[row][col] = newColor;

        // קריאה רקורסיבית לכל השכנים (למעלה, למטה, שמאלה, ימינה)
        DFS(image, row + 1, col, startColor, newColor); // למטה
        DFS(image, row - 1, col, startColor, newColor); // למעלה
        DFS(image, row, col + 1, startColor, newColor); // ימינה
        DFS(image, row, col - 1, startColor, newColor); // שמאלה
    }

    static void PrintImage(int[][] image)
    {
        foreach (var row in image)
        {
            Console.WriteLine(string.Join(" ", row));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    //חיפוש במערך דו ממדים על אי שמרומז ב 1 המציג מים ו 0 מציג יבשה על פי תנאים הבאים
    //   אי מורכב מקבוצה של תאים בעלי ערך 1 שמחוברים ביניהם אופקית או אנכית
    //  קבוצות נפרדות של תאים שמחוברים אינן נחשבות לאותו אי.

    public class NumberOfIslandsSolution
    {
        public int NumIslands(char[][] grid)
        {
            if (grid == null || grid.Length == 0) return 0;

            int numIslands = 0;

            for (int i = 0; i < grid.Length; i++)
            {
                for (int j = 0; j < grid[0].Length; j++)
                {
                    // מצאנו תא של "יבשה"
                    if (grid[i][j] == '1')
                    {
                        numIslands++;
                        DFS(grid, i, j); // הפעל חיפוש DFS כדי לסמן את כל תאי האי
                    }
                }
            }

            return numIslands;
        }

        private void DFS(char[][] grid, int row, int col)
        {
            // תנאים לבדיקת גבולות המטריצה
            if (row < 0 || col < 0 || row >= grid.Length || col >= grid[0].Length || grid[row][col] == '0')
            {
                return;
            }

            // סימון התא הנוכחי כ"ביקרנו בו"
            grid[row][col] = '0';

            // חיפוש לכל הכיוונים
            DFS(grid, row - 1, col); // למעלה
            DFS(grid, row + 1, col); // למטה
            DFS(grid, row, col - 1); // שמאלה
            DFS(grid, row, col + 1); // ימינה
        }
    }




    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    // חיפוש במטריציה על אות O בין קבוצת של X 
    // אם O מוקפת בכול הצדדים ב X נחליף אותה ל X

    public class SurroundedRegionsSolution
    {
        public void Solve(char[][] grid)
        {
            if (grid == null || grid.Length == 0 || grid[0].Length == 0) return;

            int rows = grid.Length;
            int cols = grid[0].Length;

            // פונקציית עזר לחיפוש DFS
            void DFS(int row, int col)
            {
                if (row < 0 || col < 0 || row >= rows || col >= cols || grid[row][col] != 'O')
                    return;

                grid[row][col] = 'T'; // שינוי זמני כדי לסמן תאים שאינם מוקפים
                DFS(row - 1, col); // למעלה
                DFS(row + 1, col); // למטה
                DFS(row, col - 1); // שמאלה
                DFS(row, col + 1); // ימינה
            }

            // חיפוש על גבולות המטריצה
            for (int i = 0; i < rows; i++)
            {
                DFS(i, 0); // עמודה שמאלית
                DFS(i, cols - 1); // עמודה ימנית
            }

            for (int j = 0; j < cols; j++)
            {
                DFS(0, j); // שורה עליונה
                DFS(rows - 1, j); // שורה תחתונה
            }

            // שינוי המטריצה
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (grid[i][j] == 'O')
                    {
                        grid[i][j] = 'X'; // מוקפים
                    }
                    else if (grid[i][j] == 'T')
                    {
                        grid[i][j] = 'O'; // לא מוקפים
                    }
                }
            }
        }
    }




    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    //  מבקשת למצוא את כל התאים במטריצה שמהם מים יכולים לזרום גם לאוקיינוס השקט(Pacific) וגם לאוקיינוס האטלנטי(Atlantic).



    public class PacificAtlanticSolution
    {
        public IList<IList<int>> PacificAtlantic(int[][] grid)
        {
            // רשימת תוצאות לשמירת התאים שמהם מים יכולים לזרום לשני האוקיינוסים
            var result = new List<IList<int>>();

            // בדיקת תוקף של הקלט (מטריצה ריקה או null)
            if (grid == null || grid.Length == 0 || grid[0].Length == 0)
                return result;

            // הגדרת מספר השורות והעמודות במטריצה
            int rows = grid.Length;
            int cols = grid[0].Length;

            // מטריצות בוליאניות למעקב אחר תאים שמהם מים יכולים לזרום לכל אחד מהאוקיינוסים
            bool[,] pacific = new bool[rows, cols];
            bool[,] atlantic = new bool[rows, cols];

            // פונקציית עזר לביצוע DFS
            void DFS(int row, int col, bool[,] ocean)
            {
                // אם כבר ביקרנו בתא זה, אין צורך להמשיך
                if (ocean[row, col]) return;

                // מסמנים את התא הנוכחי ככזה שמים יכולים לזרום ממנו
                ocean[row, col] = true;

                // רשימת הכיוונים האפשריים לזרימת מים: למעלה, למטה, שמאלה, ימינה
                int[][] directions = new int[][]
                {
                new int[] { -1, 0 }, // למעלה
                new int[] { 1, 0 },  // למטה
                new int[] { 0, -1 }, // שמאלה
                new int[] { 0, 1 }   // ימינה
                };

                // מעבר על כל הכיוונים
                foreach (var dir in directions)
                {
                    int newRow = row + dir[0];
                    int newCol = col + dir[1];

                    // בדיקה אם התא החדש בתוך גבולות המטריצה ואם ניתן לזרום אליו
                    if (newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols &&
                        grid[newRow][newCol] >= grid[row][col])
                    {
                        DFS(newRow, newCol, ocean); // קריאה רקורסיבית לתא החדש
                    }
                }
            }

            // ביצוע DFS עבור האוקיינוס השקט (Pacific) והאוקיינוס האטלנטי (Atlantic)

            // מעבר על השורות בגבולות השמאלי והימני
            for (int i = 0; i < rows; i++)
            {
                DFS(i, 0, pacific); // גבול שמאלי של Pacific
                DFS(i, cols - 1, atlantic); // גבול ימני של Atlantic
            }

            // מעבר על העמודות בגבולות העליון והתחתון
            for (int j = 0; j < cols; j++)
            {
                DFS(0, j, pacific); // גבול עליון של Pacific
                DFS(rows - 1, j, atlantic); // גבול תחתון של Atlantic
            }

            // איסוף תאים שמים יכולים לזרום מהם לשני האוקיינוסים
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    // אם תא מסומן בשתי המטריצות, מוסיפים אותו לתוצאות
                    if (pacific[i, j] && atlantic[i, j])
                    {
                        result.Add(new List<int> { i, j });
                    }
                }
            }

            // החזרת רשימת התוצאות
            return result;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    //סיכום

    //   עץ: DFS בעץ עובר מצומת השורש למטה עד עלה, וחוזר רק כשאין יותר ילדים להמשיך אליהם.

    void DFS(TreeNode node)
    {
        if (node == null) return;
        Console.WriteLine(node.Value); // פעולת עיבוד
        DFS(node.Left); // מעבר לתת-עץ שמאלי
        DFS(node.Right); // מעבר לתת-עץ ימני
    }

    // גרף: בגרפים יש לוודא שלא מבקרים בצומת פעמיים על ידי שימוש במערך ביקורים.


    void DFS(int node, List<int>[] graph, bool[] visited)
    {
        if (visited[node]) return;
        visited[node] = true;
        Console.WriteLine(node);
        foreach (var neighbor in graph[node])
        {
            DFS(neighbor, graph, visited);
        }
    }

    //  מטריצה: בדומה לגרף, אבל מתמקדים בשכנים בכיוונים מוגדרים מראש(למעלה, למטה, שמאלה, ימינה).
    void DFS(int row, int col, int[][] matrix, bool[,] visited, int[][] directions)
    {
        if (row < 0 || col < 0 || row >= matrix.Length || col >= matrix[0].Length || visited[row, col])
            return;
        visited[row, col] = true;
        Console.WriteLine($"Visiting: {row}, {col}");
        foreach (var dir in directions)
        {
            DFS(row + dir[0], col + dir[1], matrix, visited, directions);
        }
    }

    // 2. שימוש ב-stack ליישום איטרטיבי של DFS
    //    במקום  רקורסיה, משתמשים במחסנית כדי לעקוב אחרי הצמתים שצריך לבקר בהם.

    void DFSIterative(int[][] matrix, int startRow, int startCol)
    {
        var stack = new Stack<(int, int)>();
        var visited = new bool[matrix.Length, matrix[0].Length];
        stack.Push((startRow, startCol));

        while (stack.Count > 0)
        {
            var (row, col) = stack.Pop();
            if (row < 0 || col < 0 || row >= matrix.Length || col >= matrix[0].Length || visited[row, col])
                continue;

            visited[row, col] = true;
            Console.WriteLine($"Visiting: {row}, {col}");

            stack.Push((row - 1, col)); // Up
            stack.Push((row + 1, col)); // Down
            stack.Push((row, col - 1)); // Left
            stack.Push((row, col + 1)); // Right
        }
    }


    //  3. עבודה עם רקורסיה וניהול "עומק חיפוש"
    //  משתמשים במונה עומק כדי לעקוב אחרי השכבה הנוכחית של הרקורסיה.

    void DFSWithDepth(int node, int depth, List<int>[] graph, bool[] visited)
    {
        if (visited[node]) return;
        visited[node] = true;
        Console.WriteLine($"Node: {node}, Depth: {depth}");
        foreach (var neighbor in graph[node])
        {
            DFSWithDepth(neighbor, depth + 1, graph, visited);
        }
    }


    // מציאת מסלול בין שני צמתים.

    bool HasPath(int start, int end, List<int>[] graph, bool[] visited)
    {
        if (start == end) return true;
        if (visited[start]) return false;

        visited[start] = true;
        foreach (var neighbor in graph[start])
        {
            if (HasPath(neighbor, end, graph, visited)) return true;
        }

        return false;
    }



    //   1. פונקציית DFS גנרית


    void DFS(int node, List<int>[] graph, bool[] visited)
    {
        if (visited[node]) return;
        visited[node] = true;

        foreach (var neighbor in graph[node])
        {
            DFS(neighbor, graph, visited);
        }
    }



 //   1. מטריצות וגרפים ריקים

    if (matrix == null || matrix.Length == 0) return;














}
