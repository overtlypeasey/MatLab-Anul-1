classdef Student
    properties
        Name = [];
        ID;
        Grades = [];
    end

    methods
        function obj = Student(name, id, grades)
            obj.Name = name;
            obj.ID = id;
            obj.Grades = grades;
        end

        function obj = addGrade(obj, newGrade)
            obj.Grades = [obj.Grades newGrade];
        end

        function m = meanGrades(obj)
            m = mean(obj.Grades);
        end

        function s = getStudent(obj)
            m = meanGrades(obj);
            s = [obj.Name m];
        end
    end

end